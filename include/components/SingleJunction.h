#pragma once

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
