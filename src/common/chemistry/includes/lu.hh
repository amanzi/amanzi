/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_LU_HH_
#define AMANZI_CHEMISTRY_LU_HH_

#include <cmath>
#include <iostream>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

// prototypes
void ludcmp(double** a, int n, std::vector<int>* indx, double* d);
void lubksb(double** a, int n, const std::vector<int>& indx, std::vector<double>* b);
// void ludcmp(double** a, int n, int* indx, double* d);
// void lubksb(double** a, int n, int* indx, double b[]);

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif  // AMANZI_CHEMISTRY_LU_HH_
