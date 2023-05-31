#include <iostream>
#include "PETScAppInterface.h"
#include "DelPhiSimulation.h"

void
PETScApp::setupPETScWorkSpace()
{
  // Prepare NULL matrix
  J_Mat = NULL;
  J_MatrixFree = NULL;
  P_Mat = NULL;

  // Prepare NULL MatFDColoring
  fdcoloring = NULL;

  // Prepare PETSc vectors
  VecCreate(PETSC_COMM_WORLD, &u);
  VecSetSizes(u, PETSC_DECIDE, n_dofs);
  VecSetFromOptions(u);
  VecDuplicate(u, &u_old);
  VecDuplicate(u, &u_backup);
  VecDuplicate(u, &r);
  VecDuplicate(u, &res_tran);
  VecDuplicate(u, &res_spatial);

  // Setup SNES
  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetFunction(snes, r, snesFormFunction, (void *)(this));
  SNESMonitorSet(snes, snesMonitor, (void *)(this), NULL);

  // Setup Matrix
  setupMatrices();

  // Setup KSP/PC (most of these needs go to the inputs (with default values))
  PetscReal ksp_rtol = 1e-3;
  PetscInt ksp_maxits = 100;
  SNESGetKSP(snes, &ksp);
  KSPSetFromOptions(ksp);
  PC pc;
  KSPGetPC(ksp, &pc);
  PCFactorSetMatOrderingType(pc, MATORDERINGRCM);
  PCFactorSetLevels(pc, 5);
  KSPGMRESSetRestart(ksp, ksp_maxits);
  KSPSetTolerances(ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, ksp_maxits);

  // Finalize SNES setup
  SNESSetFromOptions(snes);
}

void
PETScApp::setupPETScIC()
{
  PetscScalar *uu;
  VecGetArray(u, &uu);
  p_sim->setupPETScIC(uu);
  VecRestoreArray(u, &uu);
}

void
PETScApp::backupSolution()
{
  VecCopy(u, u_backup);
}

void
PETScApp::restoreSolutionFromBackup()
{
  VecCopy(u_backup, u);
}

void
PETScApp::setupMatrices()
{
  // This should go to user-input
  unsigned solver_option = 2;

  if (solver_option == 0) // Newton's method + Hand-coded "(hopefully) exact" Jacobian
  {
    mooseError("Hand-calculated Jacobian is not implemented.");
    /*
    // Let the simulation setup Jacobian matrix sparsity
    p_sim->fillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc example:
    // https://petsc.org/release/src/snes/tutorials/ex1.c.html
    // Use hand-calculated Jacobian as both the Jacobian and Preconditioning Jacobian
    SNESSetJacobian(snes, P_Mat, P_Mat, formJacobian, this);*/
  }
  else if (solver_option == 1) // Matrix-free + Hand-coded Jacobian as the Preconditioning Jacobian
  {
    mooseError("Hand-calculated Jacobian is not implemented.");
    /*
    // Create Matrix-free context
    MatCreateSNESMF(snes, &J_MatrixFree);

    // Let the simulation setup Jacobian matrix sparsity
    p_sim->fillJacobianMatrixNonZeroPattern(P_Mat);

    // See PETSc example:
    // https://petsc.org/release/src/ts/tutorials/ex15.c.html
    // Use hand-calculated Jacobian as the Preconditioning Jacobian
    SNESSetJacobian(snes, J_MatrixFree, P_Mat, formJacobian, this);*/
  }
  else if (solver_option ==
           2) // Matrix-free + Finite-differencing Preconditioning Jacobian (using coloring)
  {
    // Create Matrix-free context
    MatCreateSNESMF(snes, &J_MatrixFree);

    // Let the problem setup Jacobian matrix sparsity
    p_sim->FillJacobianMatrixNonZeroEntry(P_Mat);

    // See PETSc examples:
    // https://petsc.org/release/src/snes/tutorials/ex14.c.html
    // https://petsc.org/release/src/mat/tutorials/ex16.c.html
    ISColoring iscoloring;
    MatColoring mc;
    MatColoringCreate(P_Mat, &mc);
    MatColoringSetType(mc, MATCOLORINGSL);
    MatColoringSetFromOptions(mc);
    MatColoringApply(mc, &iscoloring);
    MatColoringDestroy(&mc);
    MatFDColoringCreate(P_Mat, iscoloring, &fdcoloring);
    MatFDColoringSetFunction(
        fdcoloring, (PetscErrorCode(*)(void))(void (*)(void))snesFormFunction, this);
    MatFDColoringSetFromOptions(fdcoloring);
    MatFDColoringSetUp(P_Mat, iscoloring, fdcoloring);
    ISColoringDestroy(&iscoloring);

    SNESSetJacobian(snes,                            // snes
                    J_MatrixFree,                    // Jacobian-free
                    P_Mat,                           // Preconditioning matrix
                    SNESComputeJacobianDefaultColor, // Use finite differencing and coloring
                    fdcoloring);                     // fdcoloring
  }
  else if (solver_option == 3) // Finite-differencing, no coloring, slowest
  {
    // See PETSc example:
    // https://petsc.org/release/src/ts/tutorials/ex10.c.html
    MatCreateSeqAIJ(PETSC_COMM_SELF, n_dofs, n_dofs, PETSC_DEFAULT, PETSC_NULL, &J_Mat);
    SNESSetJacobian(snes,  // snes
                    J_Mat, // Jacobian matrix
                    J_Mat, // Preconditioning mat, use the same Jacobian mat
                    SNESComputeJacobianDefault,
                    PETSC_NULL);
  }
  else
    mooseError("Unknown Solver option.");
}

void
PETScApp::freePETScWorkSpace()
{
  // Destroy PETSc vectors
  VecDestroy(&u);
  VecDestroy(&u_old);
  VecDestroy(&r);
  VecDestroy(&res_tran);
  VecDestroy(&res_spatial);

  // Destroy PETSc matrix
  if (J_Mat != NULL)
    MatDestroy(&J_Mat);
  if (J_MatrixFree != NULL)
    MatDestroy(&J_MatrixFree);
  if (P_Mat != NULL)
    MatDestroy(&P_Mat);

  // Destroy SNES
  SNESDestroy(&snes);

  // Destroy MatFDColoring
  if (fdcoloring != NULL)
    MatFDColoringDestroy(&fdcoloring);
}

PetscErrorCode
snesFormFunction(SNES /*snes*/, Vec u, Vec f, void * appCtx)
{
  PETScApp * petscApp = (PETScApp *)appCtx;

  // zero out residuals
  VecZeroEntries(petscApp->res_tran);
  VecZeroEntries(petscApp->res_spatial);

  // get vectors
  PetscScalar *uu, *res_tran, *res_spatial;
  VecGetArray(u, &uu);
  VecGetArray(petscApp->res_tran, &res_tran);
  VecGetArray(petscApp->res_spatial, &res_spatial);

  // use the most updated solution vector to update solution, to compute RHS and transient residuals
  petscApp->p_sim->updateSolutions(uu);
  petscApp->p_sim->computeHelperVariables();
  petscApp->p_sim->computeTranRes(res_tran);
  petscApp->p_sim->computeSpatialRes(res_spatial);

  // restore vectors
  VecRestoreArray(u, &uu);
  VecRestoreArray(petscApp->res_tran, &res_tran);
  VecRestoreArray(petscApp->res_spatial, &res_spatial);

  // assemble final residuals: f = transient + spatial
  VecWAXPY(f, 1.0, petscApp->res_tran, petscApp->res_spatial);

  return 0;
}

PetscErrorCode
formJacobian(SNES /*snes*/, Vec /*u*/, Mat /*jac*/, Mat /*B*/, void * /*appCtx*/)
{
  mooseError("Hand-calculated Jacobian is not implemented.");
  /*
  PETScApp * petscApp = (PETScApp *)appCtx;
  petscApp->p_sim->computeJacobianMatrix(B);

  if (jac != B)
  {
    MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
  }*/

  return 0;
}

PetscErrorCode
snesMonitor(SNES /*snes*/, PetscInt its, PetscReal fnorm, void * /*AppCtx*/)
{
  PetscPrintf(PETSC_COMM_WORLD, "    NL Step = %2D, fnorm = %12.5E\n", its, fnorm);
  return 0;
}

PetscErrorCode
kspMonitor(KSP /*ksp*/, PetscInt its, PetscReal rnorm, void * /*AppCtx*/)
{
  PetscPrintf(PETSC_COMM_WORLD, "      Linear step = %2D, rnorm = %12.5E\n", its, rnorm);
  return 0;
}
