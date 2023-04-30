#include "CreateDelPhiSimulation.h"
#include "DelPhiMesh.h"
#include "FEProblemBase.h"
#include "CreateProblemAction.h"
#include "Factory.h"
#include "MooseApp.h"
#include "ActionWarehouse.h"

registerMooseAction("delphiApp", CreateDelPhiSimulation, "create_delphi_simulation");

InputParameters
CreateDelPhiSimulation::validParams()
{
  InputParameters params = Action::validParams();
  params.addClassDescription("Action that creates an empty DelPhiSimulation.");
  return params;
}

CreateDelPhiSimulation::CreateDelPhiSimulation(const InputParameters & params)
  : Action(params)
{
}

void
CreateDelPhiSimulation::act()
{
  // create mesh object, it is now still an empty one
  if (!_mesh)
  {
    InputParameters params = _factory.getValidParams("DelPhiMesh");
    _mesh = _factory.create<DelPhiMesh>("DelPhiMesh", "DelPhi_mesh", params);
  }

  if (!_mesh->hasMeshBase())
    _mesh->setMeshBase(_mesh->buildMeshBaseObject());

  if (!_problem)
  {
    InputParameters params = _factory.getValidParams("DelPhiSimulation");
    _app.parser().extractParams("", params); // extract global params
    // apply common parameters of the object held by CreateProblemAction to honor user inputs in
    // [Problem]
    auto p = _awh.getActionByTask<CreateProblemAction>("create_problem");
    if (p)
      params.applyParameters(p->getObjectParams());
    params.set<MooseMesh *>("mesh") = _mesh.get();
    params.set<bool>("use_nonlinear") = _app.useNonlinear();
    params.set<MooseApp *>("_moose_app") = &_app;
    _problem = _factory.create<FEProblemBase>("DelPhiSimulation", "DelPhiSimulation", params);
  }
}
