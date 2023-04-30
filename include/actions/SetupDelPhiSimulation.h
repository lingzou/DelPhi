#pragma once

#include "Action.h"

class SetupDelPhiSimulation : public Action
{
public:
  SetupDelPhiSimulation(const InputParameters & params);
  virtual ~SetupDelPhiSimulation() {}

  virtual void act() override;

public:
  static InputParameters validParams();
};
