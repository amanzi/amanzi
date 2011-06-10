/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_LU_HH_
#define AMANZI_CHEMISTRY_LU_HH_

#include <cmath>
#include <iostream>
#include <vector>

// prototypes
void ludcmp(double** a, int n, int* indx, double* d);
void lubksb(double** a, int n, int* indx, std::vector<double> &b);
void lubksb(double** a, int n, int* indx, double b[]);

#endif  // AMANZI_CHEMISTRY_LU_HH_
