#pragma once

#include "DelPhiComponent.h"

class TestComponent : public OneDComponent
{
public:
  TestComponent(const InputParameters & parameters);
  virtual ~TestComponent() {}

  virtual void buildMesh() override;
  virtual void addExternalVariables() override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

  // don't worry, TestComponent will be inherited from a different type of base component
  // for now, let's keep these functions
  virtual std::vector<OneDCell*> & getCells() override { mooseError("Error"); }
  virtual std::vector<EdgeBase*> & getEdges() override { mooseError("Error"); }
  virtual void setBoundaryEdge(DELPHI::EEndType /*end*/, EdgeBase* /*edge*/) override { mooseError("Error"); }
  virtual Real getArea() const override { mooseError("Error"); return 0.0; }

protected:
  Real _length;
  unsigned _n_elem;
  Real _dL;

  std::vector<Real> _T;
  std::vector<Real> _T_old;

  unsigned _subdomain_id;
  SubdomainName _subdomain_name;
  std::vector<Node *> _nodes;
  std::vector<Elem *> _elems;

public:
  static InputParameters validParams();
};
