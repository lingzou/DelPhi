#pragma once

#include "ExternalProblem.h"
#include "DelPhiMesh.h"
#include "PETScAppInterface.h"

class DelPhiComponent;

class DelPhiSimulation : public ExternalProblem
{
public:
  DelPhiSimulation(const InputParameters & params);
  virtual ~DelPhiSimulation();

  // Functions inherited from ExternalProblem
  virtual void externalSolve() override;
  virtual void syncSolutions(Direction direction) override;
  virtual bool converged() override;
  virtual void addExternalVariables() override;
  // Functions inherited from FEProblemBase
  virtual void onTimestepBegin() override;
  virtual void onTimestepEnd() override;

  virtual void initSimulation();
  virtual void buildMesh();

  // Residual functions - related to PETSc
  virtual void setupPETScIC(double * u);
  virtual void updateSolutions(double * u);
  virtual void computeTranRes(double * r);
  virtual void computeSpatialRes(double * r);
  // Jacobian non-zero entries - related to PETSc
  virtual void FillJacobianMatrixNonZeroEntry(Mat & P_Mat);

  virtual void
  addComponent(const std::string & type, const std::string & name, InputParameters & params);
  DelPhiMesh & phoenix_mesh() { return _delphi_mesh; }

  // Helper API functions
  void addMooseAuxVar(const std::string & name,
                      const FEType & type,
                      const std::vector<SubdomainName> & subdomain_names);

protected:
  // PETSc interface
  PETScApp * _p_PETScApp;
  // number of DOFs of this simulation
  unsigned _n_DOFs;

protected:
  struct VariableInfo
  {
    FEType _type;
    std::set<SubdomainName> _subdomain;
  };

protected:
  DelPhiMesh & _delphi_mesh;

  /// List of components in this simulation
  std::vector<MooseSharedPointer<DelPhiComponent>> _components;

  /// Variables for output purpose <var_name, sub_domain_set>
  std::map<std::string, VariableInfo> _vars;

public:
  static InputParameters validParams();
};
