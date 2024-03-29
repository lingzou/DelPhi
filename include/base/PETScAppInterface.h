#pragma once

#include <vector>
#include <set>
#include <petsc.h>

class DelPhiSimulation;

PETSC_EXTERN PetscErrorCode snesFormFunction(SNES, Vec, Vec, void *);
PETSC_EXTERN PetscErrorCode snesMonitor(SNES, PetscInt, PetscReal, void *);
PETSC_EXTERN PetscErrorCode kspMonitor(KSP, PetscInt, PetscReal, void *);
PETSC_EXTERN PetscErrorCode formJacobian(SNES, Vec, Mat, Mat, void *);

class MatrixNonZeroPattern
{
public:
  MatrixNonZeroPattern(unsigned size) { _non_zero_entries.resize(size); }
  virtual ~MatrixNonZeroPattern() {}
  virtual void addEntry(unsigned row, unsigned col) { _non_zero_entries[row].insert(col); }
  virtual void addRow(unsigned row, std::set<unsigned> & cols)
  {
    for (auto & col : cols)
      addEntry(row, col);
  }
  std::vector<std::set<unsigned>> & getNonZeroPattern() { return _non_zero_entries; }

protected:
  std::vector<std::set<unsigned>> _non_zero_entries;
};

struct PETScApp
{
  PETScApp(DelPhiSimulation * sim) : p_sim(sim) {}

  /// A pointer to the DelPhiSimulation
  DelPhiSimulation * p_sim;
  /// Number of degrees of freedoms (DOFs)
  PetscInt n_dofs;

  SNES snes;
  KSP ksp;
  MatFDColoring fdcoloring;
  /// Jacobian matrix
  Mat J_Mat;
  /// Preconditioning matrix
  Mat P_Mat;
  /// Jacobian-free matrix context
  Mat J_MatrixFree;
  /// Unknown vector
  Vec u;
  /// Unknown vector old
  Vec u_old;
  /// A backup copy for the unknown vector
  // (in case retry a failed time step, 'u' might have been into a unrecoverable status)
  Vec u_backup;
  /// Total residual (= res_transient + res_spatial)
  Vec r;
  /// Residual contribution from transient terms
  Vec res_tran;
  /// Residual contribution from non-transient spatial terms
  Vec res_spatial;

  // It is common to see a time step fail, and we have to come back to retry solving the same time step
  // with a smaller time step. Because the unknown vector 'u' might be in a unrecoverable status, we would
  // need a copy of the last known good solution.
  // The two functions are for doing backup the good solution, and get the good solution back.
  void backupSolution();
  void restoreSolutionFromBackup();

  void setupPETScWorkSpace();
  void setupPETScIC();
  void setupMatrices();
  void freePETScWorkSpace();
};
