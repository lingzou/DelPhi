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

  // data access
  virtual unsigned Ts_DOF(unsigned j, unsigned i) { return _Ts_DOFs[j][i]; }

  virtual Real Tw(unsigned j, unsigned side)
  {
    if (side == 0)        return _Tw_left[j];
    else if (side == 1)   return _Tw_right[j];
    else                  { mooseError("error"); return 0.0; }
  }

  virtual unsigned Tw_DOF(unsigned j, unsigned side)
  {
    if (side == 0)        return _Tw_left_DOFs[j]; // Tw_left_DOF_local(j) + _DOF_offset;
    else if (side == 1)   return _Tw_right_DOFs[j]; // Tw_right_DOF_local(j) + _DOF_offset;
    else                  { mooseError("error"); return 0; }
  } 

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

    // text-based output
  virtual void writeTextOutput() override;

  // don't worry, SimpleHeatStructure2 will be inherited from a different type of base component
  // for now, let's keep these functions
  virtual std::vector<OneDCell*> & getCells() override { mooseError("Error"); }
  virtual std::vector<EdgeBase*> & getEdges() override { mooseError("Error"); }
  virtual void setBoundaryEdge(DELPHI::EEndType /*end*/, EdgeBase* /*edge*/) override { mooseError("Error"); }
  virtual Real getArea() const override { mooseError("Error"); return 0.0; }

protected:
  Real _length;
  Real _width;
  Real _depth;
  unsigned _NL;
  unsigned _NW;
  Real _dL;
  Real _dW;
  Real _qv;

  const ThermalSolidProperties * _solid;

  std::vector<std::vector<Real>> _volume;
  std::vector<std::vector<Real>> _Ts;
  std::vector<std::vector<Real>> _Ts_old;
  std::vector<std::vector<Real>> _Ts_oo;

  std::vector<Real> _Tw_left;
  std::vector<Real> _Tw_right;

  std::vector<std::vector<unsigned>> _Ts_DOFs;
  std::vector<unsigned> _Tw_left_DOFs;
  std::vector<unsigned> _Tw_right_DOFs;

  // mesh related data
  unsigned _subdomain_id;
  SubdomainName _subdomain_name;
  std::vector<std::vector<Node *>> _nodes;
  std::vector<std::vector<Elem *>> _elems;

  // boundary condition related data
  Real _T_bc_left;
  Real _T_bc_right;

  OneDFlowChannel * _pipe_left;
  OneDFlowChannel * _pipe_right;
  Real _hw_left;
  Real _hw_right;

  Real _Q_internal;
  Real _Q_out_left;
  Real _Q_out_right;

  Real _rho_ref;
  Real _cp_ref;

public:
  static InputParameters validParams();
};
