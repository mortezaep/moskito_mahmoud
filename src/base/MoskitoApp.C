#include "MoskitoApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
MoskitoApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

MoskitoApp::MoskitoApp(InputParameters parameters) : MooseApp(parameters)
{
  MoskitoApp::registerAll(_factory, _action_factory, _syntax);
}

MoskitoApp::~MoskitoApp() {}

void 
MoskitoApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<MoskitoApp>(f, af, s);
  Registry::registerObjectsTo(f, {"MoskitoApp"});
  Registry::registerActionsTo(af, {"MoskitoApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
MoskitoApp::registerApps()
{
  registerApp(MoskitoApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
MoskitoApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  MoskitoApp::registerAll(f, af, s);
}
extern "C" void
MoskitoApp__registerApps()
{
  MoskitoApp::registerApps();
}
