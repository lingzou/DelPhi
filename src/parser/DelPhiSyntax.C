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
  registerTask("add_delphi_variables", true);
  registerTask("setup_delphi_ICs", true);

  try
  {
    syntax.addDependency("prepare_mesh", "create_delphi_simulation");
    syntax.addDependency("setup_mesh", "create_delphi_simulation");
    syntax.addDependency("init_delphi_simulation", "add_delphi_component");
    syntax.addDependency("build_mesh", "init_delphi_simulation");
    syntax.addDependency("init_mesh", "build_mesh");
    syntax.addDependency("add_delphi_variables", "init_mesh");
    syntax.addDependency("add_delphi_variables", "add_user_object");
    //syntax.addDependency("setup_delphi_ICs", "add_external_aux_variables");
  }
  catch (CyclicDependencyException<std::string> & e)
  {
    mooseError("Cyclic Dependency Detected during addDependency() calls");
  }
}

}
