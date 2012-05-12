/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "lu.hh"

namespace amanzi {
namespace chemistry {

static const double TINY = 1.0e-20;

void ludcmp(double** a, int n, std::vector<int>* indx, double* d) {
  int i, imax, j, k;
  double big, dum, sum, temp;
  double* vv;

  vv = new double[n];
  *d = 1.0;
  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs(a[i][j])) > big) {
        big = temp;
      }
    if (big == 0.0) {
      std::cout << "Singular matrix in routine ludcmp";
    }
    vv[i] = 1.0 / big;
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
      if ((dum = vv[i] * fabs(sum)) >= big) {
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
      *d = -*d;
      vv[imax] = vv[j];
    }
    indx->at(j) = imax;
    if (a[j][j] == 0.0) {
      a[j][j] = TINY;
    }
    if (j != n - 1) {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i < n; i++) {
        a[i][j] *= dum;
      }
    }
  }
  delete [] vv;
}  // end ludcmp()


#undef TINY

void lubksb(double** a, int n, const std::vector<int>& indx, std::vector<double>* b) {
  int i, ii = 0, ip, j;
  double sum;

  for (i = 0; i < n; i++) {
    ip = indx.at(i);
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
}  // end lubksb()

// void lubksb(double** a, int n, int* indx, double b[]) {
//   int i, ii = 0, ip, j;
//   double sum;

//   for (i = 0; i < n; i++) {
//     ip = indx[i];
//     sum = b[ip];
//     b[ip] = b[i];
//     if (ii != 0)
//       for (j = ii - 1; j < i; j++) {
//         sum -= a[i][j] * b[j];
//       }
//     else if (sum != 0.0) {
//       ii = i + 1;
//     }
//     b[i] = sum;
//   }
//   for (i = n - 1; i >= 0; i--) {
//     sum = b[i];
//     for (j = i + 1; j < n; j++) {
//       sum -= a[i][j] * b[j];
//     }
//     b[i] = sum / a[i][i];
//   }
// }  // end lubksb()

}  // namespace chemistry
}  // namespace amanzi
