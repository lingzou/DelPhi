#include "AddDelPhiComponentAction.h"
#include "DelPhiSimulation.h"

registerMooseAction("delphiApp", AddDelPhiComponentAction, "add_delphi_component");

InputParameters
AddDelPhiComponentAction::validParams()
{
  InputParameters params = MooseObjectAction::validParams();
  return params;
}

AddDelPhiComponentAction::AddDelPhiComponentAction(const InputParameters & params)
  : MooseObjectAction(params)
{
}

void
AddDelPhiComponentAction::act()
{
  DelPhiSimulation * sim_ptr = dynamic_cast<DelPhiSimulation *>(_problem.get());

  Moose::out << "Add Component: " << _type << " " << name() << std::endl;
  sim_ptr->addComponent(_type, name(), getObjectParams());
}
