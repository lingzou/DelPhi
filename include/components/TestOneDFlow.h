#pragma once

#include "DelPhiComponent.h"
#include "SinglePhaseFluidProperties.h"

class TestOneDFlow : public DelPhiComponent
{
public:
  TestOneDFlow(const InputParameters & parameters);
  virtual ~TestOneDFlow() {}

  virtual void buildMesh() override;
  virtual void addExternalVariables() override;
  virtual void setupIC(double * u) override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

protected:
  virtual void updateFluxes();

protected:
  Real _length;
  unsigned _n_elem;
  Real _dL;
  const SinglePhaseFluidProperties * _eos;
  Real _rho_ref, _rhoh_ref;

  std::vector<Real> _p;
  std::vector<Real> _p_old;
  std::vector<Real> _T;
  std::vector<Real> _T_old;
  std::vector<Real> _v;
  std::vector<Real> _v_old;

  std::vector<Real> _rho;
  std::vector<Real> _rho_old;
  std::vector<Real> _h;
  std::vector<Real> _h_old;

  std::vector<Real> _rho_edge;
  std::vector<Real> _mass_flux;
  std::vector<Real> _enthalpy_flux;

  unsigned _subdomain_id;
  SubdomainName _subdomain_name;
  std::vector<Node *> _nodes;
  std::vector<Elem *> _elems;

public:
  static InputParameters validParams();
};
