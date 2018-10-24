/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Epetra_Vector.h"

#include "DenseMatrix.hh"
#include "lapack.hh"

#include "OperatorAudit.hh"


namespace Amanzi {
namespace Operators {

/* ******************************************************************
* AAA.                                      
****************************************************************** */
int OperatorAudit::CheckSpectralBounds(int schema)
{ 
  // find location of face-based matrices
  double emin = 1e+99, emax = -1e+99;
  double cndmin = 1e+99, cndmax = 1.0, cndavg = 0.0;

  for (int i = 0; i < op_->matrices.size(); i++) {
    WhetStone::DenseMatrix Acell(op_->matrices[i]);
    int n = Acell.NumRows();

    int info, ldv(1), lwork = 4 * n;
    double VL, VR, dmem1[n], dmem2[n], dwork[lwork];
    
    WhetStone::DGEEV_F77("N", "N", &n, Acell.Values(), &n, dmem1, dmem2, 
                         &VL, &ldv, &VR, &ldv, dwork, &lwork, &info);

    OrderByIncrease_(n, dmem1);

    double e, a = dmem1[1], b = dmem1[1];  // skipping the first eigenvalue
    for (int k = 2; k < n; k++) {
      e = dmem1[k];
      a = std::min(a, e);
      b = std::max(b, e);
    }

    emin = std::min(emin, a);
    emax = std::max(emax, b);

    double cnd = fabs(b) / (fabs(a) + 1e-16);
    cndmin = std::min(cndmin, cnd);
    cndmax = std::max(cndmax, cnd);
    cndavg += cnd;
  }
  cndavg /= op_->matrices.size();

  int MyPID = 0; 
  if (MyPID == 0) {
    printf("OperatorAudit: lambda matrices\n");
    printf("   eigenvalues (min,max) = %8.3g %8.3g\n", emin, emax); 
    printf("   conditioning (min,max,avg) = %8.2g %8.2g %8.2g\n", cndmin, cndmax, cndavg);
  }
  return 0;
}


/* ******************************************************************
* Verify symmetry of the matrix.                                      
****************************************************************** */
int CheckMatrixSymmetry(Teuchos::RCP<Epetra_CrsMatrix> A)
{
  int nrows = A->NumMyRows();
  Epetra_Vector x(A->DomainMap());
  Epetra_Vector y(x), z(x);

  if (A->Comm()->getRank() == 0)
    printf("Running 10 symmetry tests: size(A)=%d\n", A->NumGlobalRows());

  for (int n = 0; n < 10; n++) {
    for (int f = 0; f < nrows; f++) {
      x[f] = double(random()) / RAND_MAX;
      y[f] = double(random()) / RAND_MAX;
    }
    double axy, ayx;
    A->Multiply(false, x, z);
    z.Dot(y, &axy);

    A->Multiply(false, y, z);
    z.Dot(x, &ayx);
    double err = fabs(axy - ayx) / (fabs(axy) + fabs(ayx) + 1e-10);
    if (A->Comm()->getRank() == 0 && err > 1e-10) {        
      printf("   Summetry violation: (Ax,y)=%12.7g (Ay,x)=%12.7g\n", axy, ayx);
    }
  }
  return 0;
}


/* ******************************************************************
* Verify coercivity of the matrix.                                      
****************************************************************** */
int CheckMatrixCoercivity(Teuchos::RCP<Epetra_CrsMatrix> A)
{
  int nrows = A->NumMyRows();
  Epetra_Vector x(A->DomainMap());
  Epetra_Vector y(x), z(x);

  if (A->Comm()->getRank() == 0)
    printf("Running 10 coercivity tests: size(A)=%d\n", A->NumGlobalRows());

  double tmp, axx, axxmin(1e+99), axxmax(-1e+99);
  for (int n = 0; n < 10; n++) {
    for (int f = 0; f < nrows; f++) {
      x[f] = double(random()) / RAND_MAX;
      y[f] = double(random()) / RAND_MAX;
    }
    A->Multiply(false, x, y);
    y.Dot(x, &axx);

    if (A->Comm()->getRank() == 0 && axx <= 1e-12) {        
      printf("   Coercivity violation: (Ax,x)=%12.7g\n", axx);
    }

    x.Norm2(&tmp);
    axx /= tmp * tmp;
    axxmin = std::min(axxmin, axx);
    axxmax = std::max(axxmax, axx);
  }

  if (A->Comm()->getRank() == 0) {        
    printf("   min/max of (Ax,x)/(x,x): %12.7g %12.7g\n", axxmin, axxmax);
  }
  return 0;
}


/* ******************************************************************
* Bubble algorithm.                                            
****************************************************************** */
void OperatorAudit::OrderByIncrease_(int n, double* mem)
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

}  // namespace Operators
}  // namespace Amanzi


