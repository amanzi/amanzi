/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water: rho = q3(T-Tref) * q1(p - pref) where
  q3 and q1 are cubic and linear polynomials, respectively.
*/

#include "H2O_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_Density::H2O_Density(Teuchos::ParameterList& eos_plist)
  : EOS_Density(eos_plist),
    ka_(999.915),
    kb_(0.0416516),
    kc_(-0.0100836),
    kd_(0.000206355),
    kT0_(273.15),
    kalpha_(5.0e-10),
    kp0_(1.0e5) {
};


double H2O_Density::Density(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_ * dT) * dT) * dT;
  return rho1bar * (1.0 + kalpha_ * (p - kp0_));
};


double H2O_Density::DDensityDT(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = kb_ + (2 * kc_ + 3 * kd_ * dT) * dT;
  return rho1bar * (1.0 + kalpha_ * (p - kp0_));
};


double H2O_Density::DDensityDp(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_ * dT) * dT) * dT;
  return rho1bar * kalpha_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi
