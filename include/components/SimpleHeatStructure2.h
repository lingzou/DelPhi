#pragma once

#include "DelPhiComponent.h"
#include "ThermalSolidProperties.h"

class SimpleHeatStructure2 : public SimpleHeatStructure
{
public:
  SimpleHeatStructure2(const InputParameters & parameters);
  virtual ~SimpleHeatStructure2() {}

  virtual void buildMesh() override;
  virtual void addExternalVariables() override;
  virtual void addPhysicalModel() override;
  virtual void setupIC(double * u) override;
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  // data access
  virtual unsigned Ts_DOF(unsigned j, unsigned i) override { return j * (_NW + 1) + i + _DOF_offset; }
  virtual Real Tw(unsigned j, unsigned side) override
  {
    if (side == 0)        return _Ts[j].front();
    else if (side == 1)   return _Ts[j].back();
    else                  { mooseError("error"); return 0.0; }
  }

  virtual unsigned Tw_DOF(unsigned j, unsigned side) override
  {
    if (side == 0)        return _Ts_DOFs[j].front();
    else if (side == 1)   return _Ts_DOFs[j].back();
    else                  { mooseError("error"); return 0; }
  }

  // Residual-related functions
  virtual void updateSolution(double * u) override;
  virtual void computeTranRes(double * r) override;
  virtual void computeSpatialRes(double * r) override;

  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp) override;

    // text-based output
  virtual void writeTextOutput() override;

  // don't worry, SimpleHeatStructure will be inherited from a different type of base component
  // for now, let's keep these functions
  virtual std::vector<OneDCell*> & getCells() override { mooseError("Error"); }
  virtual std::vector<EdgeBase*> & getEdges() override { mooseError("Error"); }
  virtual void setBoundaryEdge(DELPHI::EEndType /*end*/, EdgeBase* /*edge*/) override { mooseError("Error"); }
  virtual Real getArea() const override { mooseError("Error"); return 0.0; }

public:
  static InputParameters validParams();
};
