#pragma once

#include "GeneralPostprocessor.h"

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"
#include "OneDFlowModel.h"

class SingleJunction : public ZeroDComponent
{
public:
  SingleJunction(const InputParameters & parameters);
  virtual ~SingleJunction() {}

  virtual void addPhysicalModel() override;

protected:
  /// Name of the connecting component boundaries
  std::vector<std::string> _inputs;
  std::vector<std::string> _outputs;
  const SinglePhaseFluidProperties * _eos;

public:
  static InputParameters validParams();
};

class SNJVelocityPPS : public GeneralPostprocessor
{
public:
  SNJVelocityPPS(const InputParameters & parameters);
  virtual ~SNJVelocityPPS() {}

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual Real getValue() override { return _snj->v(); }

protected:
  snjEdge * _snj;

public:
  static InputParameters validParams();
};
