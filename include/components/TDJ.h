#pragma once

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"
#include "OneDFlowModel.h"

class TDJ : public ZeroDComponent
{
public:
  TDJ(const InputParameters & parameters);
  virtual ~TDJ() {}

  virtual void addExternalVariables() override;

protected:
  /// Name of the connecting component boundaries
  std::vector<std::string> _input;
  const SinglePhaseFluidProperties * _eos;

public:
  static InputParameters validParams();
};
