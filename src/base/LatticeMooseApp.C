#include "LatticeMooseApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
LatticeMooseApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

LatticeMooseApp::LatticeMooseApp(InputParameters parameters) : MooseApp(parameters)
{
  LatticeMooseApp::registerAll(_factory, _action_factory, _syntax);
}

LatticeMooseApp::~LatticeMooseApp() {}

void
LatticeMooseApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<LatticeMooseApp>(f, af, s);
  Registry::registerObjectsTo(f, {"LatticeMooseApp"});
  Registry::registerActionsTo(af, {"LatticeMooseApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
LatticeMooseApp::registerApps()
{
  registerApp(LatticeMooseApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
LatticeMooseApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LatticeMooseApp::registerAll(f, af, s);
}
extern "C" void
LatticeMooseApp__registerApps()
{
  LatticeMooseApp::registerApps();
}
