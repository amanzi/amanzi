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

#include "EOS_Water.hh"

namespace Amanzi {
namespace AmanziEOS {

EOS_Water::EOS_Water(Teuchos::ParameterList& eos_plist) :
    EOS_ConstantMolarMass(0.0180153),
    eos_plist_(eos_plist),

    ka_(999.915),
    kb_(0.0416516),
    kc_(-0.0100836),
    kd_(0.000206355),
    kT0_(273.15),

    kalpha_(5.0e-10),
    kp0_(1.0e5) {
};


double EOS_Water::MassDensity(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_*dT)*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_));
};


double EOS_Water::DMassDensityDT(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = kb_ + (2.0*kc_ + 3.0*kd_*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_));
};


double EOS_Water::DMassDensityDp(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_*dT)*dT)*dT;
  return rho1bar * kalpha_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi
