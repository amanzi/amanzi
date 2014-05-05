/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

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
  int m(0), nblocks = blocks_.size();

  for (int n = 0; n < nblocks; n++) {
    int schema = blocks_schema_[n];
    if (blocks_schema_[n] == schema) {
      m = n;
      break;
    }
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];

  double emin = 1e+99, emax = -1e+99;
  double cndmin = 1e+99, cndmax = 1.0, cndavg = 0.0;

  for (int i = 0; i < matrix.size(); i++) {
    WhetStone::DenseMatrix& Acell = matrix[i];
    int n = Acell.NumRows();

    int info, ldv(1), lwork = 4 * n;
    double VL, VR, dmem1[n], dmem2[n], dwork[lwork];
    
    WhetStone::DGEEV_F77("N", "N", &n, Acell.Values(), &n, dmem1, dmem2, 
                         &VL, &ldv, &VR, &ldv, dwork, &lwork, &info);

    OrderByIncrease(n, dmem1);

    double e, a = dmem1[1], b = dmem1[1];  // skipping the first eigenvalue
    for (int k = 2; k < n; k++) {
      e = dmem1[k];
      a = std::min(a, e);
      b = std::max(b, e);
    }

    emin = std::min(emin, a);
    emax = std::max(emax, b);

    double cnd = b / a;
    cndmin = std::min(cndmin, cnd);
    cndmax = std::max(cndmax, cnd);
    cndavg += cnd;
  }
  cndavg /= matrix.size();

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
int OperatorAudit::CheckMatrixSymmetry()
{
  /*
  for (int n = 0; n < 10; n++) {
    for (int f = 0; f < nfaces_owned; f++) {
      x[f] = double(random()) / RAND_MAX;
      y[f] = double(random()) / RAND_MAX;
    }
    matrix_->Aff()->Multiply(false, x, z);
    z.Dot(y, &axy);

    matrix_->Aff()->Multiply(false, y, z);
    z.Dot(x, &ayx);
    double err = fabs(axy - ayx) / (fabs(axy) + fabs(ayx) + 1e-10);
    if (MyPID == 0 && err > 1e-10) {	
      printf("   Summetry violation: (Ax,y)=%12.7g (Ay,x)=%12.7g\n", axy, ayx);
    }
  }
  */
  return 0;
}


/* ******************************************************************
* Verify coercivity of the matrix.                                      
****************************************************************** */
int OperatorAudit::CheckMatrixCoercivity()
{
  /*
  for (int n = 0; n < 10; n++) {
    for (int f = 0; f < nfaces_owned; f++) {
      x[f] = double(random()) / RAND_MAX;
      y[f] = double(random()) / RAND_MAX;
    }
    matrix_->Aff()->Multiply(false, x, y);
    y.Dot(x, &axx);

    if (MyPID == 0 && axx <= 1e-12) {	
      printf("   Coercivity violation: (Ax,x)=%12.7g\n", axx);
    }
  }
  */
  return 0;
}


/* ******************************************************************
* Bubble algorithm.                                            
****************************************************************** */
void OperatorAudit::OrderByIncrease(int n, double* mem)
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


