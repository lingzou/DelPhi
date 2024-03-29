#pragma once

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"
#include "OneDFlowModel.h"

class VolumeBranch : public ZeroDComponent
{
public:
  VolumeBranch(const InputParameters & parameters);
  virtual ~VolumeBranch();

  virtual void addPhysicalModel() override;
  virtual void setExtendedNeighbors() override;
  virtual void setupIC(double * u) override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

  // text-based output
  virtual void writeTextOutput() override;

protected:
  /// Name of the connecting component boundaries
  Real _vol;
  std::vector<std::string> _inputs;
  std::vector<std::string> _outputs;
  const SinglePhaseFluidProperties * _eos;
  Real _rho_ref, _rhoh_ref;

  BranchCell * _brh_cell;
  std::vector<brvEdgeInlet *> _outgoing_edges;
  std::vector<Real> _outgoing_areas;
  std::vector<brvEdgeOutlet *> _incoming_edges;
  std::vector<Real> _incoming_areas;

public:
  static InputParameters validParams();
};
