#pragma once

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"
#include "OneDFlowModel.h"

class TDV : public ZeroDComponent
{
public:
  TDV(const InputParameters & parameters);
  virtual ~TDV() {}

  virtual void addExternalVariables() override;

protected:
  /// Name of the connecting component boundaries
  std::vector<std::string> _input;
  const SinglePhaseFluidProperties * _eos;

public:
  static InputParameters validParams();
};
