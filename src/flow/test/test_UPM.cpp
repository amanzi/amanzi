#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "../upper_packed_matrix.hpp"
#include "../mimetic_hex.hpp"

TEST(UPM) {
  
  int n = 4;
  int m = n*(n+1)/2;

  double *a = new double[m];

  // make an tri diagonal matrix
  for (int i=0; i<m; i++) a[i] = 0.0;
  
  int l = 0;
  for (int i=0; i<n; i++) {
    l = i  +  i*(i+1)/2;
    a[l] = 2.0;
    if (l>0) a[l-1] = -1.0;
  }
  
  // store a copy 
  double *b = new double[m];
  
  for (int i=0; i<m; i++) b[i] = a[i];


  // create matrix objects, one is the inverse of the other
  upper_packed_matrix mat1 (a, n);
  mat1.invert();
  
  upper_packed_matrix mat2 (b, n);
  // now mat1 is the inverse of mat2

  double* x = new double[n];
  double* y = new double[n];
  double* aux = new double[n];

  for (int i=0; i<n; i++) {
    x[i] = 1.0;
    y[i] = 0.0;
    aux[i] = 0.0;
  }


  // check sym_matmul
  // mat2 * mat1 = indentity

  mat2.sym_matmul(x,aux);

  mat1.sym_matmul(aux,y);

  // now x and y should be the same
  CHECK_ARRAY_CLOSE(x,y,n,0.00001);

  // check col_sum
  mat2.col_sum(x);
  
  for (int i=0; i<n; i++) y[i] = 0.0;
  y[0] = 1.0; y[n-1] = 1.0;
  CHECK_ARRAY_CLOSE(x,y,n,0.00001);

  // check solve
  // since mat2 * x = aux, lets solve for x and compare 


  for (int i=0; i<n; i++) x[i] = 1.0;

  mat2.sym_matmul(x,aux);

  mat2.factor();
  mat2.solve(aux);

  CHECK_ARRAY_CLOSE(x, aux, n, 0.00001);


  delete a;
  delete b;
  delete x;
  delete y;
  delete aux;

}

