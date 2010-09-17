#ifndef __Newton_hpp__
#define __Newton_hpp__

#include "Block.hpp"

#include <iostream>
using namespace std;

class Newton {
  
public:
  Newton(int);
  virtual ~Newton();

  void solve();
  void solve(double *r, Block *J, double *conc, double *update);
  
private:

  int n;
  double *x;
  double *b;
  Block *J;
  int *indices;

};

#endif
