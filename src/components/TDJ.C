#include "TDJ.h"

registerMooseObject("delphiApp", TDJ);

InputParameters
TDJ::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addRequiredParam<Real>("v_bc", "Velocity boundary value");
  params.addRequiredParam<Real>("T_bc", "Temperature boundary value");
  params.addRequiredParam<UserObjectName>("eos", "equation of states");
  params.addRequiredParam<std::vector<std::string>>("input", "Names of the connected components");

  return params;
}

TDJ::TDJ(const InputParameters & parameters)
  : ZeroDComponent(parameters),
  _input(getParam<std::vector<std::string>>("input"))
{
  if (_input.size() != 1)
    mooseError("TDJ is expecting one connected component.");
}

void
TDJ::addExternalVariables()
{
  // handle eos first
  const UserObjectName & uo_name = getParam<UserObjectName>("eos");
  const UserObject & uo = _sim.getUserObject<UserObject>(uo_name);
  if (dynamic_cast<const SinglePhaseFluidProperties *>(&uo) == nullptr)
    mooseError("cannot convert: " + uo_name);
  else
    _eos = dynamic_cast<const SinglePhaseFluidProperties *>(&uo);

  std::string comp_name;
  DELPHI::EEndType end_type;
  DELPHI::getConnection(_input[0], comp_name, end_type);

  OneDComponent * comp_1d = dynamic_cast<OneDComponent*>(_sim.getComponentByName(comp_name));
  if (comp_1d == NULL)
    mooseError(comp_name + "is not a OneDComponent.");

  bool inlet = (end_type == DELPHI::IN);

  // boundary vEdge
  InputParameters pars = emptyInputParameters();
  pars.set<std::string>("name") = name() + ":vEdge";
  pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
  pars.set<CellBase *>("west_cell") = inlet ? NULL : (comp_1d->getCells()).back();
  pars.set<CellBase *>("east_cell") = inlet ? (comp_1d->getCells()).front():  NULL;
  pars.set<Real>("v_bc") = getParam<Real>("v_bc");
  pars.set<Real>("T_bc") = getParam<Real>("T_bc");

  EdgeBase * edge = NULL;
  if (inlet)
    edge = new vBCEdgeInlet(pars);
  else
    edge = new vBCEdgeOutlet(pars);

  comp_1d->setBoundaryEdge(end_type, edge);
}
