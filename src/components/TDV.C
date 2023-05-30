#include "TDV.h"

registerMooseObject("delphiApp", TDV);

InputParameters
TDV::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addRequiredParam<Real>("p_bc", "Pressure boundary value");
  params.addRequiredParam<Real>("T_bc", "Temperature boundary value");
  params.addRequiredParam<UserObjectName>("eos", "equation of states");
  params.addRequiredParam<std::vector<std::string>>("input", "Names of the connected components");

  return params;
}

TDV::TDV(const InputParameters & parameters)
  : ZeroDComponent(parameters),
  _input(getParam<std::vector<std::string>>("input"))
{
  if (_input.size() != 1)
    mooseError("TDV is expecting one connected component.");
}

void
TDV::addPhysicalModel()
{
  // handle eos first
  _eos = _sim.getSinglePhaseEOS(getParam<UserObjectName>("eos"));

  std::string comp_name;
  DELPHI::EEndType end_type;
  DELPHI::getConnection(_input[0], comp_name, end_type);

  OneDComponent * comp_1d = dynamic_cast<OneDComponent*>(_sim.getComponentByName(comp_name));
  if (comp_1d == NULL)
    mooseError(comp_name + "is not a OneDComponent.");

  bool inlet = (end_type == DELPHI::IN);

  // boundary pEdge
  InputParameters pars = emptyInputParameters();
  pars.set<DelPhiSimulation *>("_sim") = &_sim;
  pars.set<std::string>("name") = name() + ":pEdge";
  pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
  pars.set<CellBase *>("west_cell") = inlet ? NULL : (comp_1d->getCells()).back();
  pars.set<CellBase *>("east_cell") = inlet ? (comp_1d->getCells()).front():  NULL;
  pars.set<Real>("p_bc") = getParam<Real>("p_bc");
  pars.set<Real>("T_bc") = getParam<Real>("T_bc");

  EdgeBase * edge = NULL;
  if (inlet)
    edge = new pBCEdgeInlet(pars);
  else
    edge = new pBCEdgeOutlet(pars);

  comp_1d->setBoundaryEdge(end_type, edge);
}
