//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "MoskitoTestApp.h"
#include "MoskitoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
MoskitoTestApp::validParams()
{
  InputParameters params = MoskitoApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

MoskitoTestApp::MoskitoTestApp(InputParameters parameters) : MooseApp(parameters)
{
  MoskitoTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

MoskitoTestApp::~MoskitoTestApp() {}

void
MoskitoTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  MoskitoApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"MoskitoTestApp"});
    Registry::registerActionsTo(af, {"MoskitoTestApp"});
  }
}

void
MoskitoTestApp::registerApps()
{
  registerApp(MoskitoApp);
  registerApp(MoskitoTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
MoskitoTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MoskitoTestApp::registerAll(f, af, s);
}
extern "C" void
MoskitoTestApp__registerApps()
{
  MoskitoTestApp::registerApps();
}
