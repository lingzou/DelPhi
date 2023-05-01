#include "delphiApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"
#include "FluidPropertiesApp.h"
#include "DelPhiSyntax.h"

InputParameters
delphiApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  return params;
}

delphiApp::delphiApp(InputParameters parameters) : MooseApp(parameters)
{
  delphiApp::registerAll(_factory, _action_factory, _syntax);
}

delphiApp::~delphiApp() {}

void
delphiApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  FluidPropertiesApp::registerAll(f, af, syntax);

  Registry::registerObjectsTo(f, {"delphiApp"});
  Registry::registerActionsTo(af, {"delphiApp"});

  /* register custom execute flags, action syntax, etc. here */
  DelPhi::associateSyntax(syntax, af);
}

void
delphiApp::registerApps()
{
  registerApp(delphiApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
delphiApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  delphiApp::registerAll(f, af, s);
}
extern "C" void
delphiApp__registerApps()
{
  delphiApp::registerApps();
}
