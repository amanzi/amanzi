#ifndef __LU_hpp__
#define __LU_hpp__

#include <math.h>
#include <iostream>

// prototypes
void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

#endif 
