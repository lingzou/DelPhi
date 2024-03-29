// libmesh include
#include "libmesh/string_to_enum.h"

// moose include
#include "Executioner.h"
#include "Transient.h"

#include "DelPhiSimulation.h"
#include "DelPhiComponent.h"

registerMooseObject("delphiApp", DelPhiSimulation);

InputParameters
DelPhiSimulation::validParams()
{
  InputParameters params = ExternalProblem::validParams();
  return params;
}

DelPhiSimulation::DelPhiSimulation(const InputParameters & params)
  : ExternalProblem(params), _n_DOFs(0), _delphi_mesh(dynamic_cast<DelPhiMesh &>(_mesh))
{
  _p_PETScApp = new PETScApp(this);

  std::string fullname = getMooseApp().parser().getPrimaryFileName(false);
  size_t pos = fullname.find_last_of('.');
  fullname = fullname.substr(0, pos) + "_out.txt";
  _p_file = fopen(fullname.c_str(), "w");
}

DelPhiSimulation::~DelPhiSimulation()
{
  fclose(_p_file);
  _p_PETScApp->freePETScWorkSpace();
  delete _p_PETScApp;
}

void
DelPhiSimulation::addComponent(const std::string & type,
                                const std::string & name,
                                InputParameters & params)
{
  params.set<DelPhiSimulation *>("_sim") = this;

  MooseSharedPointer<DelPhiComponent> comp = _factory.create<DelPhiComponent>(type, name, params);
  if (_comp_by_name.find(name) == _comp_by_name.end())
    _comp_by_name[name] = comp;
  else
    mooseError("Component with name '", name, "' already exists.");

  // also put them in different 'buckets' as certain actions need to be performed in order of types of components
  if (dynamic_cast<OneDComponent*>(comp.get()) != NULL)
    _components_1d.push_back(dynamic_cast<OneDComponent*>(comp.get()));
  else if (dynamic_cast<ZeroDComponent*>(comp.get()) != NULL)
    _components_0d.push_back(dynamic_cast<ZeroDComponent*>(comp.get()));
  else
    mooseError("Can only add OneDComponent and ZeroDComponent.");
}

DelPhiComponent *
DelPhiSimulation::getComponentByName(const std::string & name)
{
  if (_comp_by_name.find(name) == _comp_by_name.end())
  {
    mooseError("Component with name '" + name + "' does not exist in the system.");
    return NULL;
  }
  else
    return (_comp_by_name.find(name))->second.get();
}

void
DelPhiSimulation::initSimulation()
{
  // before going into solving the problem, we make sure that it is a 'transient' simulation
  // using implicit-euler (BDF1) and bdf2
  Executioner * exec = _app.getExecutioner();
  if (dynamic_cast<Transient *>(exec) != NULL)
  {
    _ts = dynamic_cast<Transient *>(exec)->getTimeScheme();
    if ((_ts != Moose::TI_IMPLICIT_EULER) && (_ts != Moose::TI_BDF2))
      mooseError("Currently, only 'implicit-euler' and 'bdf2' schemes are supported.");
  }
  else
    mooseError("Currently, only 'Transiet' simulations are supported.");
}

void
DelPhiSimulation::buildMesh()
{
  for (auto & it : _comp_by_name)
    it.second->buildMesh();

  _delphi_mesh.getMesh().set_spatial_dimension(3);
  _delphi_mesh.prepare(nullptr);
}

void
DelPhiSimulation::addExternalVariables()
{
  // FEProblemBase setup PETSc options in FEProblemBase::solve()
  // From experiments, it seems to be too late to do so in ExternalProblem::solve()
  // and this is the right place to setup PETSc options
  // We may need find a better place (i.e., creating a new function) for this operation
  Moose::PetscSupport::petscSetOptions(_petsc_options, _solver_params);

  for (auto & it : _comp_by_name)
    it.second->addExternalVariables();

  // After collecting all variables info, let FEProblemBase add them
  for (auto & var : _vars)
  {
    VariableName name = var.first;
    VariableInfo & vi = var.second;
    auto order = libMesh::Utility::enum_to_string<Order>(vi._type.order);
    auto family = libMesh::Utility::enum_to_string(vi._type.family);
    std::vector<SubdomainName> subdomains(vi._subdomain.begin(), vi._subdomain.end());

    InputParameters var_params = _factory.getValidParams("MooseVariable");
    var_params.set<MooseEnum>("family") = family;
    var_params.set<MooseEnum>("order") = order;
    var_params.set<std::vector<SubdomainName>>("block") = subdomains;
    FEProblemBase::addAuxVariable("MooseVariable", name, var_params);
  }
}

void
DelPhiSimulation::addPhysicalModel()
{
  // Let each component add their variables info, there is an order here
  for (auto & comp_1d : _components_1d)
  {
    comp_1d->setDOFoffset(_n_DOFs);
    comp_1d->addPhysicalModel();
    _n_DOFs += comp_1d->getNDOF();
  }
  for (auto & comp_0d : _components_0d)
  {
    comp_0d->setDOFoffset(_n_DOFs);
    comp_0d->addPhysicalModel();
    _n_DOFs += comp_0d->getNDOF();
  }

  // Now the system is ready for preparing extended connections
  for (auto & it : _comp_by_name)
    it.second->setExtendedNeighbors();

  _p_PETScApp->n_dofs = _n_DOFs;

  // test, will figure out how to feed these info to PETSc
  // Executioner * exec = _app.getExecutioner();
  // Transient * transient = dynamic_cast<Transient *>(exec);
  // if (transient)
  //   std::cerr << "l_tol = " << transient->getParam<Real>("l_tol");

  _p_PETScApp->setupPETScWorkSpace();
  _p_PETScApp->setupPETScIC();

  Moose::out << "  Total number of nonlinear DOFs: " << _n_DOFs << std::endl;
}

const SinglePhaseFluidProperties *
DelPhiSimulation::getSinglePhaseEOS(const UserObjectName & name)
{
  const UserObject & uo = getUserObject<UserObject>(name);
  if (dynamic_cast<const SinglePhaseFluidProperties *>(&uo) == nullptr)
    mooseError("UserObject '" + name + "' is not a type of SinglePhaseFluidProperties.");
  else
    return dynamic_cast<const SinglePhaseFluidProperties *>(&uo);
}

void
DelPhiSimulation::addMooseAuxVar(const std::string & name,
                                  const FEType & type,
                                  const std::vector<SubdomainName> & subdomain_names)
{
  if (_vars.find(name) == _vars.end())
  {
    VariableInfo vi;
    vi._type = type;
    for (auto && subdomain : subdomain_names)
      vi._subdomain.insert(subdomain);

    _vars[name] = vi;
  }
  else
  {
    VariableInfo & vi = _vars[name];
    if (vi._type != type)
      mooseError(
          "A component is trying to add variable of the same name but with different order/type");
    for (auto && subdomain : subdomain_names)
      vi._subdomain.insert(subdomain);
  }

  /*
  // not working well, complaining such an IC already exist. This IC setup might be unecessary
  // since the component will write the values directly.
  {
    InputParameters params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = name;
    params.set<Real>("value") = 23.4;
    params.set<std::vector<SubdomainName>>("block") = {subdomain};
    FEProblemBase::addInitialCondition(
        "ConstantIC", name + std::string(":") + subdomain + std::string(":constIC"), params);
  }*/
}

void
DelPhiSimulation::setupPETScIC(double * u)
{
  for (auto & it : _comp_by_name)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->setupIC(u + offset);
  }
  _p_PETScApp->backupSolution();
}

void
DelPhiSimulation::computeHelperVariables()
{
  for (auto & it : _comp_by_name)
    it.second->computeHelperVariables();
}

void
DelPhiSimulation::highOrderReconstruction()
{
  for (auto & it : _comp_by_name)
    it.second->highOrderReconstruction();
}

void
DelPhiSimulation::updateSolutions(double * u)
{
  for (auto & it : _comp_by_name)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->updateSolution(u + offset);
  }
}

void
DelPhiSimulation::computeTranRes(double * r)
{
  for (auto & it : _comp_by_name)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->computeTranRes(r + offset);
  }
}

void
DelPhiSimulation::computeSpatialRes(double * r)
{
  for (auto & it : _comp_by_name)
  {
    unsigned offset = it.second->getDOFoffset();
    it.second->computeSpatialRes(r + offset);
  }
}

void
DelPhiSimulation::externalSolve()
{
  SNESSolve(_p_PETScApp->snes, NULL, _p_PETScApp->u);
}

void
DelPhiSimulation::onTimestepBegin()
{
  // onTimestepBegin() happens at the beginning of a time step before it steps into solving the time step.
  // however, sometimes a time step fails and MOOSE will try the same time step with typically a reduced time step size
  // however, onTimestepBegin() does not know if this is a fresh new time step or retry of a failed one.
  // It is safe to always get the 'good' copy of solution which stored at the end of previous successful time step
  // or initial condition if time step 0
  _p_PETScApp->restoreSolutionFromBackup();

  for (auto & it : _comp_by_name)
    it.second->onTimestepBegin();
}

void
DelPhiSimulation::onTimestepEnd()
{
  // onTimestepEnd() happens at the end of a time step and ONLY when a time step succeeds,
  // so it is a good time to keep a 'good' copy of the solution in case next time step struggles to converge and need
  // this 'good' copy to retry the time step with changes such as reduced time step size
  _p_PETScApp->backupSolution();

  // loop over components for them to handle necessary operations such as save solutions to 'old time step solutions'
  for (auto & it : _comp_by_name)
    it.second->onTimestepEnd();

  // file head for the text-based output
  fprintf(_p_file, "********************************************************************************\n");
  fprintf(_p_file, "** Time = %20.6e\n", _time);
  fprintf(_p_file, "********************************************************************************\n");
  for (auto & it : _comp_by_name)
    it.second->writeTextOutput();
}

void
DelPhiSimulation::syncSolutions(Direction /*direction*/)
{
}

bool
DelPhiSimulation::converged()
{
  // See ref: https://petsc.org/main/docs/manualpages/SNES/SNESGetConvergedReason.html
  SNESConvergedReason snes_converged_reason;
  SNESGetConvergedReason(_p_PETScApp->snes, &snes_converged_reason);
  // See ref:
  // https://petsc.org/main/docs/manualpages/SNES/SNESConvergedReason.html#SNESConvergedReason
  return (snes_converged_reason > 0);
}

void
DelPhiSimulation::FillJacobianMatrixNonZeroEntry(Mat & P_Mat)
{
  MatrixNonZeroPattern * mnzp = new MatrixNonZeroPattern(_n_DOFs);

  for (auto & it : _comp_by_name)
    it.second->FillJacobianMatrixNonZeroEntry(mnzp);

  std::vector<std::set<unsigned int>> & nzp = mnzp->getNonZeroPattern();
  PetscInt * nnz;
  PetscMalloc(_n_DOFs * sizeof(PetscInt), &nnz);
  for (unsigned i = 0; i < nzp.size(); i++)
    nnz[i] = nzp[i].size();

  MatCreateSeqAIJ(PETSC_COMM_SELF, _n_DOFs, _n_DOFs, 0, nnz, &P_Mat);

  PetscReal zero = 0.0;
  for (unsigned i = 0; i < nzp.size(); i++) // loop on rows
  {
    PetscInt row = i;
    for (auto j : nzp[i]) // for each row, loop on its columns
    {
      PetscInt col = j;
      MatSetValues(P_Mat, 1, &row, 1, &col, &zero, INSERT_VALUES);
    }
  }
  MatAssemblyBegin(P_Mat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P_Mat, MAT_FINAL_ASSEMBLY);

  // MatView(P_Mat, PETSC_VIEWER_STDOUT_SELF);
  // MatView(P_Mat, PETSC_VIEWER_DRAW_WORLD);

  delete mnzp;
  PetscFree(nnz);
}
