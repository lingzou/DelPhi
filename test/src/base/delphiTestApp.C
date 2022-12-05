//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "delphiTestApp.h"
#include "delphiApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
delphiTestApp::validParams()
{
  InputParameters params = delphiApp::validParams();
  return params;
}

delphiTestApp::delphiTestApp(InputParameters parameters) : MooseApp(parameters)
{
  delphiTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

delphiTestApp::~delphiTestApp() {}

void
delphiTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  delphiApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"delphiTestApp"});
    Registry::registerActionsTo(af, {"delphiTestApp"});
  }
}

void
delphiTestApp::registerApps()
{
  registerApp(delphiApp);
  registerApp(delphiTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
delphiTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  delphiTestApp::registerAll(f, af, s);
}
extern "C" void
delphiTestApp__registerApps()
{
  delphiTestApp::registerApps();
}
