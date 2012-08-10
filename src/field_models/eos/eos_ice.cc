/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid ice.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_ice.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldModel,EOSIce> EOSIce::factory_("EOS of ice");

EOSIce::EOSIce(Teuchos::ParameterList& eos_plist) :
    EOS(eos_plist),

    ka_(916.724),
    kb_(-0.147143),
    kc_(-0.000238095),

    kT0_(273.15),

    kalpha_(1.0e-10),
    kp0_(1.0e5) {

  InitializeFromPlist_();
};


EOSIce::EOSIce(const EOSIce& other) :
    EOS(other),
    ka_(916.724),
    kb_(-0.147143),
    kc_(-0.000238095),
    kT0_(273.15),
    kalpha_(1.0e-10),
    kp0_(1.0e5),
    M_(other.M_) {}

// ---------------------------------------------------------------------------
// Virtual copy constructor.
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldModel> EOSIce::Clone() const {
  return Teuchos::rcp(new EOSIce(*this));
}

double EOSIce::Density(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + kc_*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_)) / M_;
};

double EOSIce::DDensityDT(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = kb_ + 2.0*kc_*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_)) / M_;
};

double EOSIce::DDensityDp(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + kc_*dT)*dT;
  return rho1bar * kalpha_ / M_;
};


void EOSIce::InitializeFromPlist_() {
  if (eos_plist_.isParameter("Molar mass of ice [kg/mol]")) {
    M_ = eos_plist_.get<double>("Molar mass of ice [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("Molar mass of ice [g/mol]", 18.0153)*1e-3;
  }
};

} // namespace
} // namespace
} // namespace
