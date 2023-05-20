#include "VolumeBranch.h"

registerMooseObject("delphiApp", VolumeBranch);

InputParameters
VolumeBranch::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addParam<Real>("volume", 0.0, "Volume of the volume branch");
  params.addRequiredParam<UserObjectName>("eos", "equation of states");
  params.addRequiredParam<std::vector<std::string>>("inputs", "Names of the inlet component");
  params.addRequiredParam<std::vector<std::string>>("outputs", "Names of the outlet component");

  params.addRequiredParam<Real>("initial_P", "Initial value for pressure");
  params.addRequiredParam<Real>("initial_T", "Initial value for temperature");

  return params;
}

VolumeBranch::VolumeBranch(const InputParameters & parameters)
  : ZeroDComponent(parameters),
  _vol(getParam<Real>("volume")),
  _inputs(getParam<std::vector<std::string>>("inputs")),
  _outputs(getParam<std::vector<std::string>>("outputs"))
{
}

VolumeBranch::~VolumeBranch()
{
  if (_brh_cell) delete _brh_cell;
}

void
VolumeBranch::setExtendedNeighbors()
{
  _brh_cell->setExtendedNeighborCells();
}

void
VolumeBranch::addExternalVariables()
{
  _n_DOFs = 2; // p and T of this branch

  // handle eos first
  _eos = _sim.getSinglePhaseEOS(getParam<UserObjectName>("eos"));

  // BranchCell of this VolumeBranch
  InputParameters pars = emptyInputParameters();
  pars.set<std::string>("name") = name() + ":BranchCell";
  pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
  pars.set<Real>("dL") = 0; // we will revisit this value

  _brh_cell = new BranchCell(pars);
  _brh_cell->setDOF(_DOF_offset, _DOF_offset + 1);

  // inputs or outputs do not really matter, so join them together for easier looping
  std::vector<std::string> connections;
  connections.insert(connections.end(), _inputs.begin(), _inputs.end());
  connections.insert(connections.end(), _outputs.begin(), _outputs.end());

  for (unsigned i = 0; i < connections.size(); i++)
  {
    std::string comp_name;
    DELPHI::EEndType end_type;
    DELPHI::getConnection(connections[i], comp_name, end_type);

    OneDComponent * comp_1d = dynamic_cast<OneDComponent*>(_sim.getComponentByName(comp_name));
    if (comp_1d == NULL)
      mooseError(comp_name + "is not a OneDComponent.");

    if (end_type == DELPHI::IN)
    {
      // brvEdgeInlet
      InputParameters pars = emptyInputParameters();
      pars.set<std::string>("name") = name() + ":brvEdgeInlet:" + std::to_string(i);
      pars.set<CellBase *>("west_cell") = _brh_cell;
      pars.set<CellBase *>("east_cell") = (comp_1d->getCells()).front();
      pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;

      brvEdgeInlet * edge = new brvEdgeInlet(pars);
      //_brh_cell->addEdge(edge, 1.0);
      _outgoing_edges.push_back(edge);
      _outgoing_areas.push_back(comp_1d->getArea());
      comp_1d->setBoundaryEdge(end_type, edge);
    }
    else if (end_type == DELPHI::OUT)
    {
      // brvEdgeOutlet
      InputParameters pars = emptyInputParameters();
      pars.set<std::string>("name") = name() + ":brvEdgeOutlet:" + std::to_string(i);
      pars.set<CellBase *>("west_cell") = (comp_1d->getCells()).back();
      pars.set<CellBase *>("east_cell") =  _brh_cell;
      pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;

      brvEdgeOutlet * edge = new brvEdgeOutlet(pars);
      //_brh_cell->addEdge(edge, -1.0);
      _incoming_edges.push_back(edge);
      _incoming_areas.push_back(comp_1d->getArea());
      comp_1d->setBoundaryEdge(end_type, edge);
    }
    else
      mooseError("Unknown connection type.");
  }
}

void
VolumeBranch::setupIC(double * u)
{
  Real p_init = getParam<Real>("initial_P");
  Real T_init = getParam<Real>("initial_T");

  u[0] = p_init;
  u[1] = T_init;
  _brh_cell->initialize(p_init, T_init);

  _rho_ref = _eos->rho_from_p_T(p_init, T_init);
  _rhoh_ref = _rho_ref * _eos->h_from_p_T(p_init, T_init);
}

void
VolumeBranch::updateSolution(double * u)
{
  _brh_cell->updateSolution(u[0], u[1]);
}

void
VolumeBranch::computeTranRes(double * res)
{
  res[0] = _vol * (_brh_cell->rho() - _brh_cell->rho_o()) / _sim.dt() / _rho_ref;
  res[1] = _vol * (_brh_cell->rhoh() - _brh_cell->rhoh_o()) / _sim.dt() / _rhoh_ref;
}

void
VolumeBranch::computeSpatialRes(double * res)
{
  for (unsigned i = 0; i < _outgoing_edges.size(); i++)
  {
    Real area = _outgoing_areas[i];

    res[0] += _outgoing_edges[i]->mass_flux() * area; // rho * u * A
    res[1] += _outgoing_edges[i]->enthalpy_flux() * area; // rho * u * h * A
  }
  for (unsigned i = 0; i < _incoming_edges.size(); i++)
  {
    Real area = _incoming_areas[i];

    res[0] -= _incoming_edges[i]->mass_flux() * area;
    res[1] -= _incoming_edges[i]->enthalpy_flux() * area;
  }
  // scaling
  res[0] /= _rho_ref;
  res[1] /= _rhoh_ref;
}

void
VolumeBranch::onTimestepBegin()
{
}

void
VolumeBranch::onTimestepEnd()
{
  // save old solutions
  _brh_cell->saveOldSlns();
}

void
VolumeBranch::writeTextOutput()
{
  FILE * file = _sim.getTextOutputFile();

  fprintf(file, "Component = %s\n", name().c_str());
  fprintf(file, "%20s%20s%20s%20s%20s\n", "Cell", "p", "T", "rho", "h");
  fprintf(file, "%20s%20.8e%20.8e%20.8e%20.8e\n", _brh_cell->name().c_str(), _brh_cell->p(), _brh_cell->T(), _brh_cell->rho(), _brh_cell->h());
}

void
VolumeBranch::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  mnzp->addRow(_brh_cell->pDOF(), _brh_cell->getConnectedDOFs());
  mnzp->addRow(_brh_cell->TDOF(), _brh_cell->getConnectedDOFs());
}
