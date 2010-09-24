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
  
 private:

  int n;
  double *x;
  double *r;
  Block *J;
  int *indices;

};

#endif // __Newton_hpp__
