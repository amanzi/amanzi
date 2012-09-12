/*
This is the audit component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_SerialDenseVector.h"

#include "Matrix_MFD.hpp"
#include "Matrix_Audit.hpp"

namespace Amanzi {

/* ******************************************************************
* Constructor.                                      
****************************************************************** */
Matrix_Audit::Matrix_Audit(Teuchos::RCP<AmanziMesh::Mesh> mesh, AmanziFlow::Matrix_MFD* matrix)
{ 
  matrix_type = MATRIX_AUDIT_MFD; 
  mesh_ = mesh;
  matrix_ = matrix;
}


/* ******************************************************************
* Destructor.                                      
****************************************************************** */
Matrix_Audit::~Matrix_Audit()
{
  delete [] dmem1;
  delete [] dmem2;
  delete [] dwork1;
}


/* ******************************************************************
* Calculates global information about matrices                                       
****************************************************************** */
void Matrix_Audit::InitAudit()
{
  printf("Matrix_Audit: initializing for matrix id =%2d\n", matrix_type); 
  MyPID = mesh_->cell_map(false).Comm().MyPID();

  if (matrix_type == MATRIX_AUDIT_MFD) {
    A = &(matrix_->Aff_cells());
  }

  lda = 1;
  for (int i = 0; i < A->size(); i++) {
    Teuchos::SerialDenseMatrix<int, double>& Ai = (*A)[i];
    lda = std::max<int>(lda, Ai.numRows());
  }

  // allocate memory for Lapack
  dmem1 = new double[lda + 1];
  dmem2 = new double[lda + 1];

  lwork1 = 10 * (lda + 1);
  dwork1 = new double[lwork1];
  printf("Matrix_Audit: maximum matrix size =%3d\n", lda); 
}


/* ******************************************************************
* AAA.                                      
****************************************************************** */
int Matrix_Audit::RunAudit()
{
  int ierr;
  ierr = CheckSpectralBounds();
  ierr |= CheckSpectralBoundsExtended();
  ierr |= CheckSpectralBoundsSchurComplement();
  return ierr;
}


/* ******************************************************************
* AAA.                                      
****************************************************************** */
int Matrix_Audit::CheckSpectralBounds()
{ 
  Teuchos::LAPACK<int, double> lapack;
  int info;
  double VL, VR;

  double emin = 1e+99, emax = -1e+99;
  double cndmin = 1e+99, cndmax = 1.0, cndavg = 0.0;

  for (int c = 0; c < A->size(); c++) {
    Teuchos::SerialDenseMatrix<int, double> Acell((*A)[c]);
    int n = Acell.numRows();
    
    lapack.GEEV('N', 'N', n, Acell.values(), n, dmem1, dmem2, 
                &VL, 1, &VR, 1, dwork1, lwork1, &info);

    OrderByIncrease(n, dmem1);

    double e, a = dmem1[1], b = dmem1[1];  // skipping the first eigenvalue
    for (int k=2; k<n; k++) {
      e = dmem1[k];
      if (e > 0.99) continue;
      a = std::min<double>(a, e);
      b = std::max<double>(b, e);
    }

    emin = std::min<double>(emin, a);
    emax = std::max<double>(emax, b);

    double cnd = b / a;
    cndmin = std::min<double>(cndmin, cnd);
    cndmax = std::max<double>(cndmax, cnd);
    cndavg += cnd;
  }
  cndavg /= A->size();

  if (MyPID == 0) {
    printf("Matrix_Audit: lambda matrices\n");
    printf("   eigenvalues (min,max) = %8.3g %8.3g\n", emin, emax); 
    printf("   conditioning (min,max,avg) = %8.2g %8.2g %8.2g\n", cndmin, cndmax, cndavg);
  }
  return 0;
}


/* ******************************************************************
* AAA.                                      
****************************************************************** */
int Matrix_Audit::CheckSpectralBoundsExtended()
{ 
  Errors::Message msg;

  Teuchos::LAPACK<int, double> lapack;
  int info;
  double VL, VR;

  double emin = 1e+99, emax = -1e+99;
  double cndmin = 1e+99, cndmax = 1.0, cndavg = 0.0;

  for (int c = 0; c < A->size(); c++) {
    Teuchos::SerialDenseMatrix<int, double>& Aff = matrix_->Aff_cells()[c];
    Epetra_SerialDenseVector& Acf = matrix_->Acf_cells()[c];
    Epetra_SerialDenseVector& Afc = matrix_->Afc_cells()[c];
    double& Acc = matrix_->Acc_cells()[c];
    int n = Aff.numRows();

    Teuchos::SerialDenseMatrix<int, double> Acell(n+1, n+1);
    for (int i = 0; i < n; i++) {
      Acell(n, n) = Acc;
      Acell(i, n) = Afc[i];
      Acell(n, i) = Acf[i];
      for (int j = 0; j < n; j++) Acell(i, j) = Aff(i, j);
    }

    if (Acc <= 0.0) {
      cout << Acell << endl;
      msg << "Matrix Audit: Acc is not positive.";
      Exceptions::amanzi_throw(msg);
    }
    
    n++;
    lapack.GEEV('N', 'N', n, Acell.values(), n, dmem1, dmem2, 
                &VL, 1, &VR, 1, dwork1, lwork1, &info);

    OrderByIncrease(n, dmem1);

    double e, a = dmem1[1], b = dmem1[1];  // skipping the first eigenvalue
    for (int k=2; k<n; k++) {
      e = dmem1[k];
      if (e > 0.99) continue;
      a = std::min<double>(a, e);
      b = std::max<double>(b, e);
    }

    emin = std::min<double>(emin, a);
    emax = std::max<double>(emax, b);

    double cnd = b / a;
    cndmin = std::min<double>(cndmin, cnd);
    cndmax = std::max<double>(cndmax, cnd);
    cndavg += cnd;
  }
  cndavg /= A->size();

  if (MyPID == 0) {
    printf("Matrix_Audit: p-lambda matrices\n");
    printf("   eigenvalues (min,max) = %8.3g %8.3g\n", emin, emax); 
    printf("   conditioning (min,max,avg) = %8.2g %8.2g %8.2g\n", cndmin, cndmax, cndavg);
  }
  return 0;
}


/* ******************************************************************
* AAA.                                      
****************************************************************** */
int Matrix_Audit::CheckSpectralBoundsSchurComplement()
{ 
  Teuchos::LAPACK<int, double> lapack;
  int info;
  double VL, VR;

  double emin = 1e+99, emax = -1e+99;
  double cndmin = 1e+99, cndmax = 1.0, cndavg = 0.0;

  for (int c = 0; c < A->size(); c++) {
    Teuchos::SerialDenseMatrix<int, double> Acell((*A)[c]);
    Epetra_SerialDenseVector& Acf = matrix_->Acf_cells()[c];
    Epetra_SerialDenseVector& Afc = matrix_->Afc_cells()[c];
    double& Acc = matrix_->Acc_cells()[c];
    int n = Acell.numRows();

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Acell(i, j) -= Afc[i] * Acf[j] / Acc;
      }
    }

    lapack.GEEV('N', 'N', n, Acell.values(), n, dmem1, dmem2, 
                &VL, 1, &VR, 1, dwork1, lwork1, &info);

    OrderByIncrease(n, dmem1);

    double e, a = dmem1[1], b = dmem1[1];  // skipping the first eigenvalue
    for (int k=2; k<n; k++) {
      e = dmem1[k];
      if (e > 0.99) continue;
      a = std::min<double>(a, e);
      b = std::max<double>(b, e);
    }

    emin = std::min<double>(emin, a);
    emax = std::max<double>(emax, b);

    double cnd = b / a;
    cndmin = std::min<double>(cndmin, cnd);
    cndmax = std::max<double>(cndmax, cnd);
    cndavg += cnd;
  }
  cndavg /= A->size();

  if (MyPID == 0) {
    printf("Matrix_Audit: Schur complement matrices\n");
    printf("   eigenvalues (min,max) = %8.3g %8.3g\n", emin, emax); 
    printf("   conditioning (min,max,avg) = %8.2g %8.2g %8.2g\n", cndmin, cndmax, cndavg);
  }
  return 0;
}


/* ******************************************************************
* Bubble algorithm.                                            
****************************************************************** */
void Matrix_Audit::OrderByIncrease(int n, double* mem)
{
  for (int i = 0; i < n; i++) {
    for (int j = 1; j < n-i; j++) {
      if (mem[j-1] > mem[j]) {
         double tmp = mem[j];
         mem[j] = mem[j-1];
         mem[j-1] = tmp;
      }
    }
  }
}

}  // namespace Amanzi


