#include "LU.hpp"
#include "Newton.hpp"

Newton::Newton(const int n_) {
  n = n_;
  J = new Block(n);
  x = new double[n_];
  b = new double[n_];
  indices = new int[n_];
}

void Newton::solve() {
  cout << "Solved!\n";
}

void Newton::solve(double *r, Block *J, double *conc, double *update) {
  
  // scale the Jacobian
  for (int i=0; i<n; i++) {
    double max = J->getRowAbsMax(i);
    if (max > 1.) {
      double scale = 1./max;
      b[i] = r[i]*scale;
      J->scaleRow(i,scale);
    }
  }

  // for derivatives with respect to ln concentration
  for (int i=0; i<n; i++)
    J->scaleRow(i,conc[i]);

  double D;
  ludcmp(J->getValues(),n,indices,&D);
  lubksb(J->getValues(),n,indices,b);
}

Newton::~Newton() {
  if (J) delete J;
  if (x) delete [] x;
  if (b) delete [] b;
  if (indices) delete [] indices;
}
