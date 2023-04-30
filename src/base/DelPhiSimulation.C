#include "libmesh/string_to_enum.h"

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
}

DelPhiSimulation::~DelPhiSimulation()
{
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
  _components.push_back(comp);
}

void
DelPhiSimulation::initSimulation()
{
}

void
DelPhiSimulation::buildMesh()
{
  for (auto & comp : _components)
    comp->buildMesh();

  _delphi_mesh.getMesh().set_spatial_dimension(3);
  _delphi_mesh.prepare(true);
}

void
DelPhiSimulation::addExternalVariables()
{
  // Let each component add their variables info
  for (auto & comp : _components)
  {
    comp->setDOFoffset(_n_DOFs);
    comp->addExternalVariables();
    _n_DOFs += comp->getNDOF();
  }
  _p_PETScApp->n_dofs = _n_DOFs;
  _p_PETScApp->setupPETScWorkSpace();

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
DelPhiSimulation::updateSolutions(double * u)
{
  for (auto & comp : _components)
  {
    unsigned offset = comp->getDOFoffset();
    comp->updateSolution(u + offset);
  }
}

void
DelPhiSimulation::computeTranRes(double * r)
{
  for (auto & comp : _components)
  {
    unsigned offset = comp->getDOFoffset();
    comp->computeTranRes(r + offset);
  }
}

void
DelPhiSimulation::computeSpatialRes(double * r)
{
  for (auto & comp : _components)
  {
    unsigned offset = comp->getDOFoffset();
    comp->computeSpatialRes(r + offset);
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
  for (auto & comp : _components)
    comp->onTimestepBegin();
}

void
DelPhiSimulation::onTimestepEnd()
{
  for (auto & comp : _components)
    comp->onTimestepEnd();
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

  for (auto & comp : _components)
    comp->FillJacobianMatrixNonZeroEntry(mnzp);

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
