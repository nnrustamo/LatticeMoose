//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "LatticeMooseTestApp.h"
#include "LatticeMooseApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
LatticeMooseTestApp::validParams()
{
  InputParameters params = LatticeMooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

LatticeMooseTestApp::LatticeMooseTestApp(InputParameters parameters) : MooseApp(parameters)
{
  LatticeMooseTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

LatticeMooseTestApp::~LatticeMooseTestApp() {}

void
LatticeMooseTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  LatticeMooseApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"LatticeMooseTestApp"});
    Registry::registerActionsTo(af, {"LatticeMooseTestApp"});
  }
}

void
LatticeMooseTestApp::registerApps()
{
  registerApp(LatticeMooseApp);
  registerApp(LatticeMooseTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
LatticeMooseTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LatticeMooseTestApp::registerAll(f, af, s);
}
extern "C" void
LatticeMooseTestApp__registerApps()
{
  LatticeMooseTestApp::registerApps();
}
