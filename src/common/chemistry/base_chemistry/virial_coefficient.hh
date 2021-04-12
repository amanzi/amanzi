/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#ifndef AMANZI_CHEMISTRY_VIRIAL_COEFFICIENT_HH_
#define AMANZI_CHEMISTRY_VIRIAL_COEFFICIENT_HH_

#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

class VirialCoefficient {
 public:
  VirialCoefficient();
  ~VirialCoefficient() {};

  void UpdateVirial(const double& temp, const double& pressure);

  double GetVirial() const { return virial; }

  void SetPol(double poli) {
    npol++;
    pol.push_back(poli);
  }

  int GetIsp1() const { return isp1; }
  int GetIsp2() const { return isp2; }
  int GetIsp3() const { return isp3; }

  int GetIfun1() const { return ifun1; }
  int GetIfun2() const { return ifun2; }
  int GetIfun3() const { return ifun3; }

  void SetIsp1(const int isp1_) { isp1 = isp1_; }
  void SetIsp2(const int isp2_) { isp2 = isp2_; }
  void SetIsp3(const int isp3_) { isp3 = isp3_; }

  void SetIfun1(const int ifun1_) { ifun1 = ifun1_; }
  void SetIfun2(const int ifun2_) { ifun2 = ifun2_; }
  void SetIfun3(const int ifun3_) { ifun3 = ifun3_; }

 private:
  std::vector<double> pol;

  int npol;
  double virial;

  int isp1, isp2, isp3;
  int ifun1, ifun2, ifun3;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
