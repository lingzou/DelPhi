#pragma once

#include "MooseObjectAction.h"

/**
 * The action creates DelPhi Components and includes them in the Simlation.
 */
class AddDelPhiComponentAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  AddDelPhiComponentAction(const InputParameters & params);
  virtual ~AddDelPhiComponentAction() {}

  virtual void act() override;
};
