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
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactoryWithState<ViscosityModel,ViscosityModelWater> ViscosityModelWater::factory_("liquid water");

ViscosityModelWater::ViscosityModelWater(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S) :
    ViscosityModel(eos_plist, S),
    kav1_(998.333),
    kbv1_(-8.1855),
    kcv1_(0.00585),
    kbv2_(1.3272),
    kcv2_(-0.001053),
    kT1_(293.15) {

  InitializeFromPlist_();
};


ViscosityModelWater::ViscosityModelWater(const ViscosityModelWater& other) :
    ViscosityModel(other),
    kav1_(998.333),
    kbv1_(-8.1855),
    kcv1_(0.00585),
    kbv2_(1.3272),
    kcv2_(-0.001053),
    kT1_(293.15) {}

// ---------------------------------------------------------------------------
// Virtual copy constructor.
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldModel> ViscosityModelWater::Clone() const {
  return Teuchos::rcp(new ViscosityModelWater(*this));
}


double ViscosityModelWater::Viscosity(double T) {
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

double ViscosityModelWater::DViscosityDT(double T) {
  double dT = kT1_ - T;
  double xi;

  ASSERT(0);
};



void ViscosityModelWater::InitializeFromPlist_() {};

} // namespace
} // namespace
} // namespace
