#pragma once

#include "DelPhiComponent.h"
#include "ThermalSolidProperties.h"

class SimpleHeatStructure : public OneDComponent
{
public:
  SimpleHeatStructure(const InputParameters & parameters);
  virtual ~SimpleHeatStructure() {}

  virtual void buildMesh() override;
  virtual void addExternalVariables() override;
  virtual void addPhysicalModel() override;
  virtual void setupIC(double * u) override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

  // don't worry, SimpleHeatStructure will be inherited from a different type of base component
  // for now, let's keep these functions
  virtual std::vector<OneDCell*> & getCells() override { mooseError("Error"); }
  virtual std::vector<EdgeBase*> & getEdges() override { mooseError("Error"); }
  virtual void setBoundaryEdge(DELPHI::EEndType /*end*/, EdgeBase* /*edge*/) override { mooseError("Error"); }
  virtual Real getArea() const override { mooseError("Error"); return 0.0; }

protected:
  Real _length;
  Real _width;
  unsigned _NL;
  unsigned _NW;
  Real _dL;
  Real _dW;

  const ThermalSolidProperties * _solid;

  std::vector<std::vector<Real>> _volume;
  std::vector<std::vector<Real>> _Ts;
  std::vector<std::vector<Real>> _Ts_old;
  std::vector<std::vector<Real>> _Ts_oo;

  unsigned _subdomain_id;
  SubdomainName _subdomain_name;
  std::vector<std::vector<Node *>> _nodes;
  std::vector<std::vector<Elem *>> _elems;

public:
  static InputParameters validParams();
};
