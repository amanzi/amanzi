/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_water.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<EOS,EOSWater> EOSWater::factory_("liquid water");

EOSWater::EOSWater(Teuchos::ParameterList& eos_plist) :
    eos_plist_(eos_plist),

    ka_(999.915),
    kb_(0.0416516),
    kc_(0.0100836),
    kd_(0.000206355),
    kT0_(273.15),

    kalpha_(5.0e-10),
    kp0_(1.0e5),

    kav1_(998.333),
    kbv1_(-8.1855),
    kcv1_(0.00585),
    kbv2_(1.3272),
    kcv2_(-0.001053),
    kT1_(293.15) {

  InitializeFromPlist_();
};

double EOSWater::MassDensity(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_*dT)*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_));
};

double EOSWater::DMassDensityDT(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = kb_ + (2.0*kc_ + 3.0*kd_*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_));

};

double EOSWater::DMassDensityDp(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_*dT)*dT)*dT;
  return rho1bar * kalpha_;
};

double EOSWater::MolarDensity(double T, double p) {
  return MassDensity(T,p) / M_;
};

double EOSWater::DMolarDensityDT(double T, double p) {
  return DMassDensityDT(T,p) / M_;
};

double EOSWater::DMolarDensityDp(double T, double p) {
  return DMassDensityDp(T,p) / M_;
};

double EOSWater::Viscosity(double T) {
  double dT = kT1_ - T;
  double xi;
  if (T < kT1_) {
    double A = kav1_ + (kbv1_ + kcv1_*dT)*dT;
    xi = 1301.0 * (1.0/A - 1.0/kav1_);
  } else {
    double A = (kbv2_ + kcv2_*dT)*dT;
    xi = A/(T - 168.15);
  }
  return 0.001 * pow(10.0, xi);
};

void EOSWater::InitializeFromPlist_() {
  if (eos_plist_.isParameter("Molar mass of water [kg/mol]")) {
    M_ = eos_plist_.get<double>("Molar mass of water [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("Molar mass of water [g/mol]", 18.0153)*1e-3;
  }
};

} // namespace
} // namespace
} // namespace
