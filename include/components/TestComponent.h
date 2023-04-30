#pragma once

#include "DelPhiComponent.h"

class TestComponent : public DelPhiComponent
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
