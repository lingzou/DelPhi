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
    _mesh(_sim.phoenix_mesh()),
    _n_DOFs(0),
    _DOF_offset(0)
{
}
