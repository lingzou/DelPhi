#pragma once

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"
#include "OneDFlowModel.h"

class TestOneDFlow : public OneDComponent
{
public:
  TestOneDFlow(const InputParameters & parameters);
  virtual ~TestOneDFlow();

  virtual void buildMesh() override;
  virtual void addExternalVariables() override;
  virtual void setExtendedNeighbors() override;
  virtual void setupIC(double * u) override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

  // data access
  virtual std::vector<OneDCell*> & getCells() override { return _cells; }
  virtual std::vector<EdgeBase*> & getEdges() override { return _edges; }
  virtual void setBoundaryEdge(DELPHI::EEndType end, EdgeBase* edge) override;
  virtual Real getArea() const override { return _flow_area; }

  // text-based output
  virtual void writeTextOutput() override;

protected:
  Real _length;
  unsigned _n_elem;
  Real _dL;
  Real _flow_area;
  Real _dh;
  Real _f;
  Real _qv;
  const SinglePhaseFluidProperties * _eos;
  Real _rho_ref, _rhoh_ref;

  std::vector<OneDCell*> _cells;
  std::vector<EdgeBase*> _edges;

  unsigned _subdomain_id;
  SubdomainName _subdomain_name;
  std::vector<Node *> _nodes;
  std::vector<Elem *> _elems;

public:
  static InputParameters validParams();
};
