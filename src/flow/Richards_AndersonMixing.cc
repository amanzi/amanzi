/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Author:  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Local routine
****************************************************************** */
void AndersonMixingMatrix(
    Teuchos::SerialDenseMatrix<int, double>& A, Epetra_MultiVector& krylov, 
    Epetra_MultiVector& d_krylov, int index, int m) 
{
  if (m == 0) return;

  // populate matrix
  double dot;
  Epetra_Vector* v1 = d_krylov(index);
  for (int i = 0; i < m; i++) {
    Epetra_Vector* v2 = d_krylov(i);
    v1->Dot(*v2, &dot);
    A(i, index) = A(index, i) = dot; 
    A(i, m) = 1.0;
    A(m, i) = -1.0;
  }

  // create right-hand side
  Teuchos::SerialDenseVector<int, double> b(m + 1);
  for (int i = 0; i < m; i++) {
    A(m, i) = 1.0;
    A(i, m) = -1.0;
    b(i) = 0.0;
  }
  b(m) = 1.0;

  // solve the problem
  Teuchos::SerialDenseMatrix<int, double> Acopy(A);

  Teuchos::LAPACK<int, double> lapack;
  int info, lda = Acopy.numRows();
  int ipiv[m + 1];

  lapack.GESV(m + 1, 1, Acopy.values(), lda, ipiv, b.values(), m + 1, &info);
  if (info != 0) {
    Errors::Message msg;
    msg << "Flow PK: Anderson acceleration failed in Lapack.";
    Exceptions::amanzi_throw(msg);
  }

  // update the 
  v1 = krylov(index);
  for (int n = 0; n < v1->MyLength(); n++) {
    double sum = 0.0;
    for( int i = 0; i < m; i++) {
      Epetra_Vector* v2 = krylov(i);
      sum += (*v2)[n] * b(i);
    }
    (*v1)[n] = sum;
  }
}


/* ******************************************************************
* Makes one AA step during transient time integration.
* This is the experimental method.     
****************************************************************** */
int Richards_PK::AndersonMixingTimeStep(double Tp, double dTp, double& dTnext)
{
  int mmax = 3;  // maximum number of Krylov vectors
  // allocate memory
  Epetra_MultiVector krylov(solution->Map(), mmax);
  Epetra_MultiVector d_krylov(solution->Map(), mmax);

  Teuchos::SerialDenseMatrix<int, double> A(mmax + 1, mmax + 1);

  // create solver
  if (is_matrix_symmetric) 
      solver->SetAztecOption(AZ_solver, AZ_cg);
  else
      solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, AZ_none);
  solver->SetAztecOption(AZ_conv, AZ_rhs);

  // initialize solver
  Epetra_Vector* solution_new = krylov(mmax - 1);
  Epetra_Vector* solution_old = krylov(0);
  *solution_new = *solution;
  *solution_old = *solution;

  int itrs, m = 0, index = mmax - 1;
  for (itrs = 0; itrs < 20; itrs++) {
    // update the last Krylov vector
    AndersonMixingMatrix(A, krylov, d_krylov, index, m);

    // update pointers
    solution_old = solution_new;
    index = (++index) % mmax;
    solution_new = krylov(index);

    Epetra_Vector* solution_old_cells = FS->CreateCellView(*solution_old);
    Epetra_Vector* solution_new_faces = FS->CreateFaceView(*solution_new);

    // create algebraic problem
    CalculateRelativePermeability(*solution_new);

    double time = Tp + dTp;
    UpdateSourceBoundaryData(time, *solution_new_faces);

    matrix_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
    matrix_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, matrix_);
    AddTimeDerivative_MFDpicard(*solution, *solution_cells, dTp, matrix_);
    matrix_->ApplyBoundaryConditions(bc_model, bc_values);
    matrix_->AssembleGlobalMatrices();
    rhs = matrix_->rhs();

    // create preconditioner
    preconditioner_->CreateMFDstiffnessMatrices(*Krel_cells, *Krel_faces, Krel_method);
    preconditioner_->CreateMFDrhsVectors();
    AddGravityFluxes_MFD(K, *Krel_cells, *Krel_faces, Krel_method, preconditioner_);
    AddTimeDerivative_MFDpicard(*solution, *solution_cells, dTp, preconditioner_);
    preconditioner_->ApplyBoundaryConditions(bc_model, bc_values);
    preconditioner_->AssembleGlobalMatrices();
    preconditioner_->ComputeSchurComplement(bc_model, bc_values);
    preconditioner_->UpdatePreconditioner();

    // call AztecOO solver
    solver->SetRHS(&*rhs);
    solver->SetLHS(&*solution_new);

    solver->Iterate(max_itrs_linear, convergence_tol_linear);
    int num_itrs = solver->NumIters();

    // update d_krylov history
    m = std::min<int>(++m, mmax);
    Epetra_Vector* solution_diff = d_krylov(index);
    *solution_diff = *solution_new;
    solution_diff->Update(-1.0, *solution_old, 1.0);

    // error analysis
    Epetra_Vector* solution_diff_cells = FS->CreateCellView(*solution_diff);
    double error = ErrorNormPicardExperimental(*solution_old_cells, *solution_diff_cells);

    if (MyPID == 0) { // && verbosity >= FLOW_VERBOSITY_HIGH) {
      double linear_residual = solver->ScaledResidual();
      std::printf("Picard:%4d   Pressure(error=%9.4e)  solver(%8.3e, %4d)\n",
          itrs, error, linear_residual, num_itrs);
    }

    if (error < 1e-4 && itrs > 0) break;
  }

  // make decision on the next time step
  solution_new = krylov(index);
  if (itrs < 10) {
    *solution = *solution_new;
    dTnext = min(1e+4, dT * 1.2);
  } else if (itrs < 15) {
    *solution = *solution_new;
    dTnext = dT;
  } else if (itrs < 19) {
    *solution = *solution_new;
    dTnext = dT * 0.5;
    throw itrs;
  }

  return itrs;
}


}  // namespace AmanziFlow
}  // namespace Amanzi

