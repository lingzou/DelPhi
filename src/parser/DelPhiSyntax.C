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

  registerTask("init_delphi_simulation", true);
  registerTask("build_mesh", true);
  registerTask("add_physical_model", true);
  registerTask("setup_delphi_ICs", true);

  try
  {
    syntax.addDependency("prepare_mesh", "create_delphi_simulation");
    syntax.addDependency("setup_mesh", "create_delphi_simulation");
    syntax.addDependency("build_mesh", "add_delphi_component");
    syntax.addDependency("init_mesh", "build_mesh");
    syntax.addDependency("add_physical_model", "add_function");
  }
  catch (CyclicDependencyException<std::string> & e)
  {
    mooseError("Cyclic Dependency Detected during addDependency() calls");
  }
}

}
