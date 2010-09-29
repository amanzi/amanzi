// adapted from the F90 module 
// UPPER_PACKED_MATRIX that was written
// by Neil N. Carlson on 10 Jan 2006
//
// adapted by M. Berndt 14 Sep 2010


#include "upper_packed_matrix.hpp"
#include "math.h"

void upper_packed_matrix::factor()
{

  int l, m, p, q;
  double s, t;

  // Compute the upper Choleski factor.

  a[0] = sqrt(a[0]);
  p = 1;
  while (p < n*(n+1)/2 ) {
    l = 0;
    q = p;
    t = 0.0;
    while (l < p) {
      s = a[q];
      m = p;
      while (m < q) {
	s = s - a[l]*a[m];
	l = l + 1;
	m = m + 1;
      }
      s = s / a[l];
      a[q] = s;
      t = t + s*s;
      l = l + 1;
      q = q + 1;
      if (q >= n*(n+1)/2 ) s = sqrt(double( n*(n+1)/2 - q - 1));  // trigger an exception
    }
    s = a[q] - t;
    a[q] = sqrt(s);
    p = q + 1;
  }
}


void upper_packed_matrix::invert() {

  int l, m, p, q;
  double s, t;
  
  if (n == 0) return;

  // Compute the upper Choleski factor.
  factor();

  // Compute the inverse of the Choleski factor.
  a[0] = 1.0 / a[0];
  p = 1;
  while (p < n*(n+1)/2 ) {
    l = 0;
    q = p;
    while (l < p) {
      s = a[q];
      m = p;
      while (m < q) {
	a[m] = a[m] + s*a[l];
	l = l + 1;
	m = m + 1;
      }
      a[m] = s*a[l];
      l = l + 1;
      q = q + 1;
    }
    s = 1.0 / a[q];
    a[q] = s;
    s = -s;
    m = p;
    while (m < q) {
      a[m] = s*a[m];
      m = m + 1;
    }
    p = q + 1;
  }

  // Compute the product of the inverse Choleski factors.
  a[0] = a[0]*a[0];
  p = 1;
  while (p < n*(n+1)/2 ) {
    l = 0;
    q = p;
    while (l < p) {
      s = a[q];
      m = p;
      while (m <= q) {
	a[l] = a[l] + s*a[m];
	l = l + 1;
	m = m + 1;
      }
      q = q + 1;
    }
    s = a[q];
    while (l <= q) {
      a[l] = s*a[l];
      l = l + 1;
    }
    p = q + 1;
  }

}


void upper_packed_matrix::solve(double *x) {
  
  int i, j, l;
  double s;
  
  
  // Forward substitution.
  x[0] = x[0]/a[0];
  l = 1;
  for (i = 1; i<n; i++) {
    s = x[i];
    for (j = 0; j<=i-1; j++) { 
      s = s - a[l]*x[j];
      l = l + 1;
    }
    x[i] = s/a[l];
    l = l + 1;
  }
  
  // Backward substitution.
  l = l - 1;
  for (i = n-1; i>=1; i--) {
    x[i] = x[i]/a[l];
    l = l - 1;
    for (j = i-1; j>=0; j--) {
      x[j] = x[j] - a[l]*x[i];
      l = l - 1;
    }
  }
  x[0] = x[0]/a[0];
}


void upper_packed_matrix::col_sum(double *csum) {
  
  int i, j, l;
  double s;
    
  l = 0;
  for (i = 0; i<n; i++) {
    s = 0.0;
    for (j = 0; j<i; j++) {
      s = s + a[l];
      csum[j] = csum[j] + a[l];
      l = l + 1;
    }
    csum[i] = s + a[l];
    l = l + 1;
  }
   
}


void upper_packed_matrix::sym_matmul(double *b, double *c) {

  int i, j, l;
  double bi, ci;
  
  l = 0;
  for (i = 0; i<n; i++) {
    bi = b[i];
    ci = 0.0;
    for (j = 0; j<=i-1; j++) {
      ci = ci + a[l]*b[j];
      c[j] = c[j] + a[l] * bi;
      l = l + 1;
    }
    c[i] = ci + a[l]*bi;
    l = l + 1;
  }

}
