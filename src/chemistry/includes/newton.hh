#ifndef AMANZI_CHEMISTRY_NEWTON_HH_
#define AMANZI_CHEMISTRY_NEWTON_HH_

#include "block.hh"

#include <iostream>
#include <vector>
#include <cmath>

class Newton {
  
public:
  Newton(int);
  virtual ~Newton();

  void LUDecomposition(double **a, int n, int *indx);
  void LUBackSolve(double **a, int n, int *indx, std::vector<double> &b);

  void size(int i) { this->size_ = i; }
  int size(void) const { return this->size_; }

  void solve();
  
 private:

  int size_;
  std::vector<double> x_;
  std::vector<double> r_;
  Block *J_;

  double d_;
  std::vector<int> indices_;
  std::vector<double> vv_;

};

#endif // AMANZI_CHEMISTRY_NEWTON_HH_
