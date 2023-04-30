#pragma once

#include "Action.h"

class CreateDelPhiSimulation : public Action
{
public:
  CreateDelPhiSimulation(const InputParameters & params);
  virtual ~CreateDelPhiSimulation() {}

  virtual void act() override;

public:
  static InputParameters validParams();
};
