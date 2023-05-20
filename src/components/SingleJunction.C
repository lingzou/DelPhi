#include "SingleJunction.h"

registerMooseObject("delphiApp", SingleJunction);

InputParameters
SingleJunction::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addRequiredParam<UserObjectName>("eos", "equation of states");
  params.addRequiredParam<std::vector<std::string>>("inputs", "Names of the inlet component");
  params.addRequiredParam<std::vector<std::string>>("outputs", "Names of the outlet component");

  return params;
}

SingleJunction::SingleJunction(const InputParameters & parameters)
  : ZeroDComponent(parameters),
  _inputs(getParam<std::vector<std::string>>("inputs")),
  _outputs(getParam<std::vector<std::string>>("outputs"))
{
  if (_inputs.size() != 1)
    mooseError("SingleJunction is expecting one inlet component.");
  if (_outputs.size() != 1)
    mooseError("SingleJunction is expecting one outlet component.");
}

void
SingleJunction::addExternalVariables()
{
  // handle eos first
  _eos = _sim.getSinglePhaseEOS(getParam<UserObjectName>("eos"));

  // component at 'input' connection
  std::string comp_name_in;
  DELPHI::EEndType end_type_in;
  DELPHI::getConnection(_inputs[0], comp_name_in, end_type_in);

  OneDComponent * comp_1d_in = dynamic_cast<OneDComponent*>(_sim.getComponentByName(comp_name_in));
  if (comp_1d_in == NULL)
    mooseError(comp_name_in + "is not a OneDComponent.");

  // component at 'output' connection
  std::string comp_name_out;
  DELPHI::EEndType end_type_out;
  DELPHI::getConnection(_outputs[0], comp_name_out, end_type_out);

  OneDComponent * comp_1d_out = dynamic_cast<OneDComponent*>(_sim.getComponentByName(comp_name_out));
  if (comp_1d_out == NULL)
    mooseError(comp_name_out + "is not a OneDComponent.");

  // let's make sure that SingleJunction is connecting the outlet (out) of one pipe and
  //   the inlet (in) of the other pipe (order does not matter)
  if (end_type_in == end_type_out)
    mooseError("For SingleJunction, '" + name() + "', please connect it with the (out) of one pipe"
                " and the (in) of the other pipe, e.g.,\n  inputs = 'pipe-1(out)'\n  outputs = 'pipe-2(in)'\n"
                "The order does not matter, e.g.,\n  outputs = 'pipe-1(out)'\n  inputs = 'pipe-2(in)'\n");

  // snjEdge
  InputParameters pars = emptyInputParameters();
  pars.set<std::string>("name") = name() + ":snjEdge";
  pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
  if (end_type_in == DELPHI::OUT) // pipe_in (out) <-> SNJ <-> pipe_out (in)
  {
    pars.set<CellBase *>("west_cell") = (comp_1d_in->getCells()).back();
    pars.set<CellBase *>("east_cell") = (comp_1d_out->getCells()).front();
  }
  else // pipe_in (in) <-> SNJ <-> pipe_out (out)
  {
    pars.set<CellBase *>("west_cell") = (comp_1d_out->getCells()).back();
    pars.set<CellBase *>("east_cell") = (comp_1d_in->getCells()).front();
  }

  snjEdge * edge = new snjEdge(pars);
  comp_1d_in->setBoundaryEdge(end_type_in, edge);

  InputParameters pars_s = emptyInputParameters();
  pars_s.set<std::string>("name") = name() + ":snjEdge_shadow";
  pars_s.set<snjEdge *>("real_edge") = edge;
  pars_s.set<const SinglePhaseFluidProperties *>("eos") = _eos;
  // a show edge just mimic what the true edge does but having no real connections and v_bc, T_bc
  pars_s.set<CellBase *>("west_cell") = NULL;
  pars_s.set<CellBase *>("east_cell") = NULL;
  pars_s.set<Real>("v_bc") = 0.0; // dummy value
  pars_s.set<Real>("T_bc") = 0.0; // dummy value
  EdgeBase * shadow_edge = new snjShadowEdge(pars_s);
  comp_1d_out->setBoundaryEdge(end_type_out, shadow_edge);
}
