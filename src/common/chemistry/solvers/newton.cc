/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include "newton.hh"

namespace Amanzi {
namespace AmanziChemistry {

Newton::Newton(const int n) {
  size(n);
  x_.resize(n);
  r_.resize(n);
  indices_.resize(n);
  vv_.resize(n);
}

void Newton::solve() {
  std::cout << "Solved!\n";
}


void Newton::LUDecomposition(double** a, int n, int* indx) {
  const double small_number = 1.e-20;
  int i, imax, j, k;
  double big, dum, sum, temp;

  d_ = 1.0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = std::fabs(a[i][j])) > big) {
        big = temp;
      }
    if (big == 0.0) {
      std::cout << "Singular matrix in routine ludcmp";
    }
    vv_[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++) {
    for (i = 0; i < j; i++) {
      sum = a[i][j];
      for (k = 0; k < i; k++) {
        sum -= a[i][k] * a[k][j];
      }
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++) {
      sum = a[i][j];
      for (k = 0; k < j; k++) {
        sum -= a[i][k] * a[k][j];
      }
      a[i][j] = sum;
      if ((dum = vv_[i] * std::fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 0; k < n; k++) {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      d_ = -d_;
      vv_[imax] = vv_[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0) {
      a[j][j] = small_number;
    }
    if (j != n - 1) {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i < n; i++) {
        a[i][j] *= dum;
      }
    }
  }
}

#undef TINY

void Newton::LUBackSolve(double** a, int n, int* indx, std::vector<double>* b) {
  int i, ii = 0, ip, j;
  double sum;

  for (i = 0; i < n; i++) {
    ip = indx[i];
    sum = b->at(ip);
    (*b)[ip] = b->at(i);
    if (ii != 0) {
      for (j = ii - 1; j < i; j++) {
        sum -= a[i][j] * b->at(j);
      }
    } else if (sum != 0.0) {
      ii = i + 1;
    }
    (*b)[i] = sum;
  }
  for (i = n - 1; i >= 0; i--) {
    sum = b->at(i);
    for (j = i + 1; j < n; j++) {
      sum -= a[i][j] * b->at(j);
    }
    (*b)[i] = sum / a[i][i];
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
