#include "SetupDelPhiSimulation.h"

registerMooseAction("delphiApp", SetupDelPhiSimulation, "init_delphi_simulation");
registerMooseAction("delphiApp", SetupDelPhiSimulation, "build_mesh");
registerMooseAction("delphiApp", SetupDelPhiSimulation, "add_physical_model");
//registerMooseAction("delphiApp", SetupDelPhiSimulation, "setup_delphi_ICs");

InputParameters
SetupDelPhiSimulation::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Actions associated with setting up DelPhiSimulation.");
  return params;
}

SetupDelPhiSimulation::SetupDelPhiSimulation(const InputParameters & params) : Action(params) {}

void
SetupDelPhiSimulation::act()
{
  DelPhiSimulation * sim_ptr = dynamic_cast<DelPhiSimulation *>(_problem.get());
  if (sim_ptr == NULL) // This is not a DelPhi simulation, e.g., running a moose input file
    return;

  if (_current_task == "init_delphi_simulation")
  {
    Moose::out << "Initialize DelPhi Simulation:" << std::endl;
    sim_ptr->initSimulation();
  }
  else if (_current_task == "build_mesh")
  {
    Moose::out << "Create mesh... ";
    sim_ptr->buildMesh();
    Moose::out << "OK" << std::endl;
  }
  else if (_current_task == "add_physical_model")
  {
    Moose::out << "Add physical model... ";
    sim_ptr->addPhysicalModel();
    Moose::out << "OK" << std::endl;
  }
  /*
  else if (_current_task == "setup_delphi_ICs")
  {
    Moose::out << "Setup PETSc initial conditions... ";
    sim_ptr->setupIC();
    Moose::out << "OK" << std::endl;
  }*/
  else
    mooseError("Unknown task associated with SetupDelPhiSimulation: " + _current_task);
}
