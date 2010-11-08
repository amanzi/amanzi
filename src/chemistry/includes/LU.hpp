#ifndef __LU_hpp__
#define __LU_hpp__

#include <cmath>
#include <iostream>
#include <vector>

// prototypes
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, std::vector<double> &b);
void lubksb(double **a, int n, int *indx, double b[]);

#endif // __LU_hpp__
