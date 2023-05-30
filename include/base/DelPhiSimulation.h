#pragma once

#include "ExternalProblem.h"
#include "DelPhiMesh.h"
#include "PETScAppInterface.h"
#include "SinglePhaseFluidProperties.h"

class DelPhiComponent;
class OneDComponent;
class ZeroDComponent;

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
  virtual void addPhysicalModel();

  // Residual functions - related to PETSc
  virtual void setupPETScIC(double * u);
  virtual void updateSolutions(double * u);
  virtual void computeTranRes(double * r);
  virtual void computeSpatialRes(double * r);
  // Jacobian non-zero entries - related to PETSc
  virtual void FillJacobianMatrixNonZeroEntry(Mat & P_Mat);

  virtual void
  addComponent(const std::string & type, const std::string & name, InputParameters & params);
  virtual DelPhiComponent * getComponentByName(const std::string & name);
  virtual DelPhiMesh & phoenix_mesh() { return _delphi_mesh; }

  // Helper API functions
  virtual const SinglePhaseFluidProperties * getSinglePhaseEOS(const UserObjectName & name);
  virtual void addMooseAuxVar(const std::string & name,
                              const FEType & type,
                              const std::vector<SubdomainName> & subdomain_names);
  // will work on MOOSE-stype output file later
  virtual FILE * getTextOutputFile() { return _p_file; }

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

  /// Map of components by their names, [name, comp_ptr]
  std::map<std::string, MooseSharedPointer<DelPhiComponent>> _comp_by_name;
  /// components in different 'buckets'
  std::vector<OneDComponent*> _components_1d;
  std::vector<ZeroDComponent*> _components_0d;

  /// Variables for output purpose <var_name, sub_domain_set>
  std::map<std::string, VariableInfo> _vars;

  /// text-based output
  FILE * _p_file;

public:
  static InputParameters validParams();
};
