#include "DelPhiSyntax.h"
#include "ActionFactory.h"
#include "Syntax.h"

namespace DelPhi
{

void
associateSyntax(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  syntax.registerActionSyntax("CreateDelPhiSimulation", "Components", "create_delphi_simulation");
  registerTask("create_delphi_simulation", false);

  registerSyntaxTask("AddDelPhiComponentAction", "Components/*", "add_delphi_component");
  registerMooseObjectTask("add_delphi_component", Component, false);

  // With or without the [GlobalSimParameters] input block, 'setup_global_sim_parameters' action will be executed
  //   Without the [GlobalSimParameters] input block, or an empty [GlobalSimParameters] block
  //   Default global simulation parameters will be used
  syntax.registerActionSyntax("SetupGlobalSimParameters", "GlobalSimParameters", "setup_global_sim_parameters");
  registerTask("setup_global_sim_parameters", true);

  registerTask("init_delphi_simulation", true);
  registerTask("build_mesh", true);
  registerTask("add_physical_model", true);
  registerTask("setup_delphi_ICs", true);

  try
  {
    // clang-format off
    // order of execution: top to bottom
    syntax.addDependencySets("(create_delphi_simulation)"
                             "(setup_global_sim_parameters)"
                             "(setup_mesh)"
                             "(add_delphi_component)"
                             "(build_mesh)"
                             "(init_mesh)"
                             "(setup_executioner)"
                             "(add_physical_model)");
    // clang-format on
  }
  catch (CyclicDependencyException<std::string> & e)
  {
    mooseError("Cyclic Dependency Detected during addDependency() calls");
  }
}

}
