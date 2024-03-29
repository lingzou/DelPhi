#include "DelPhiComponent.h"

unsigned DelPhiComponent::subdomain_id = 0;

InputParameters
DelPhiComponent::validParams()
{
  InputParameters params = MooseObject::validParams();
  params.addPrivateParam<DelPhiSimulation *>("_sim");
  params.addPrivateParam<std::string>("built_by_action", "add_phoenix_component");
  params.registerBase("Component");

  return params;
}

DelPhiComponent::DelPhiComponent(const InputParameters & parameters)
  : MooseObject(parameters),
    _sim(*getParam<DelPhiSimulation *>("_sim")),
    _mesh(_sim.delphi_mesh()),
    _n_DOFs(0),
    _DOF_offset(0)
{
}

InputParameters
OneDComponent::validParams()
{
  InputParameters params = DelPhiComponent::validParams();
  return params;
}

InputParameters
ZeroDComponent::validParams()
{
  InputParameters params = DelPhiComponent::validParams();
  return params;
}
