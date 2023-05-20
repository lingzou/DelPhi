#pragma once

#include "MooseObject.h"
#include "DelPhiSimulation.h"
#include "DelPhiMesh.h"
#include "DelPhiTypes.h"
#include "OneDFlowModel.h"

class DelPhiComponent : public MooseObject
{
public:
  DelPhiComponent(const InputParameters & parameters);
  virtual ~DelPhiComponent() {}

  virtual void buildMesh()
  { /*not all components have mesh*/
  }
  virtual void addExternalVariables()
  { /*not all components have variables*/
  }
  virtual void setExtendedNeighbors()
  {
  }
  virtual void setupIC(double * /*u*/)
  { /*not all components have initial conditions*/
  }
  virtual void onTimestepBegin()
  { /*not all components need act at time step begin*/
  }
  virtual void onTimestepEnd()
  { /*not all components need act at time step end*/
  }

  virtual void setDOFoffset(unsigned offset) final { _DOF_offset = offset; }
  virtual unsigned getNDOF() const final { return _n_DOFs; }
  virtual unsigned getDOFoffset() const final { return _DOF_offset; }

  // Residual-related functions
  virtual void updateSolution(double * /*u*/)
  { /*not all components have variables*/
  }
  virtual void computeTranRes(double * /*r*/)
  { /*not all components have residuals*/
  }
  virtual void computeSpatialRes(double * /*r*/)
  { /*not all components have residuals*/
  }
  virtual void FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * /*mnzp*/)
  { /*not all components have residuals thus Jacobian*/
  }

  virtual void writeTextOutput() {}

protected:
  DelPhiSimulation & _sim;
  DelPhiMesh & _mesh;

protected:
  // a counter to get the current subdomain_id
  static unsigned subdomain_id;
  unsigned getNextSubdomainId() { return subdomain_id++; }

  // number of DOF this component owns
  unsigned _n_DOFs;
  // the first DOF in the global simulation
  unsigned _DOF_offset;

public:
  static InputParameters validParams();
};

class OneDComponent : public DelPhiComponent
{
public:
  OneDComponent(const InputParameters & parameters) : DelPhiComponent(parameters) {}
  virtual ~OneDComponent() {}

  // data access
  virtual std::vector<OneDCell*> & getCells() = 0;
  virtual std::vector<EdgeBase*> & getEdges() = 0;
  virtual void setBoundaryEdge(DELPHI::EEndType end, EdgeBase* edge) = 0;
  virtual Real getArea() const = 0;

public:
  static InputParameters validParams();
};

class ZeroDComponent : public DelPhiComponent
{
public:
  ZeroDComponent(const InputParameters & parameters) : DelPhiComponent(parameters) {}
  virtual ~ZeroDComponent() {}

public:
  static InputParameters validParams();
};
