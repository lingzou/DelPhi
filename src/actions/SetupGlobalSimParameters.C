#include "SetupGlobalSimParameters.h"

registerMooseAction("delphiApp", SetupGlobalSimParameters, "setup_global_sim_parameters");

InputParameters
SetupGlobalSimParameters::validParams()
{
  InputParameters params = Action::validParams();

  // gravity
  VectorValue<Real> default_gravity = VectorValue<Real>(0.0, 0.0, -9.81);
  params.addParam<RealVectorValue>("gravity", default_gravity, "gravity");

  params.addClassDescription("Actions to set up global DelPhiSimulation input parameters.");
  return params;
}

SetupGlobalSimParameters::SetupGlobalSimParameters(const InputParameters & params)
  : Action(params)
{
}

void
SetupGlobalSimParameters::act()
{
  DelPhiSimulation * sim = dynamic_cast<DelPhiSimulation*>(_problem.get());
  if (sim)
  {
    InputParameters & pars = const_cast<InputParameters &>(sim->parameters());
    // All GlobalSimParameters are sent to DelPhiSimulation's parameter list
    pars += _pars;
  }
}
