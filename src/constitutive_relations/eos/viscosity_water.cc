/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "viscosity_water.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<ViscosityRelation,ViscosityWater> ViscosityWater::factory_("liquid water");

ViscosityWater::ViscosityWater(Teuchos::ParameterList& eos_plist) :
    eos_plist_(eos_plist),
    kav1_(998.333),
    kbv1_(-8.1855),
    kcv1_(0.00585),
    kbv2_(1.3272),
    kcv2_(-0.001053),
    kT1_(293.15) {};


double ViscosityWater::Viscosity(double T) {
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

double ViscosityWater::DViscosityDT(double T) {
  double dT = kT1_ - T;
  double xi;

  ASSERT(0);
};


} // namespace
} // namespace
