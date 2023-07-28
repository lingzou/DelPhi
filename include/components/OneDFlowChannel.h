#pragma once

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"
#include "OneDFlowModel.h"

class SimpleHeatStructure;

class OneDFlowChannel : public OneDComponent
{
public:
  OneDFlowChannel(const InputParameters & parameters);
  virtual ~OneDFlowChannel();

  virtual void buildMesh() override;
  virtual void addExternalVariables() override;
  virtual void addPhysicalModel() override;
  virtual void setExtendedNeighbors() override;
  virtual void setupIC(double * u) override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void highOrderReconstruction() override;
  virtual void computeHelperVariables() override;
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

  // conjugate heat transfer
  virtual void addWallHeating(SimpleHeatStructure * hs, Real hw, Real aw, unsigned side)
  {
    // unsigned side will be enum
    _hs.push_back(hs); _hs_side.push_back(side);
    _hws.push_back(hw); _aws.push_back(aw);
  }

protected:
  unsigned _order;

  Real _length;
  unsigned _n_elem;
  Real _dL;
  Real _gL;
  Real _flow_area;
  Real _dh;
  Real _qv;
  bool _has_Tw;
  Real _Tw;
  Real _hw;
  Real _aw;
  const SinglePhaseFluidProperties * _eos;
  Real _rho_ref, _rhoh_ref;

  std::vector<OneDCell*> _cells;
  std::vector<EdgeBase*> _edges;

  unsigned _subdomain_id;
  SubdomainName _subdomain_name;
  std::vector<Node *> _nodes;
  std::vector<Elem *> _elems;

  std::vector<SimpleHeatStructure *> _hs;
  std::vector<unsigned> _hs_side;
  std::vector<Real> _hws;
  std::vector<Real> _aws;

public:
  static InputParameters validParams();
};
