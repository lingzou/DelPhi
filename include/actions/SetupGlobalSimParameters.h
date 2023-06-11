#pragma once

#include "Action.h"

class SetupGlobalSimParameters : public Action
{
public:
  SetupGlobalSimParameters(const InputParameters & params);
  virtual ~SetupGlobalSimParameters() {}

  virtual void act() override;

public:
  static InputParameters validParams();
};
