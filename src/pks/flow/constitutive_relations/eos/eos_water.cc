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
Utils::RegisteredFactoryWithState<EOS,EOSWater> EOSWater::factory_("liquid water");

EOSWater::EOSWater(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S) :
    EOS(eos_plist, S),

    ka_(999.915),
    kb_(0.0416516),
    kc_(0.0100836),
    kd_(0.000206355),
    kT0_(273.15),

    kalpha_(5.0e-10),
    kp0_(1.0e5) {
  InitializeFromPlist_();
};


EOSWater::EOSWater(const EOSWater& other) :
    EOS(other),
    ka_(999.915),
    kb_(0.0416516),
    kc_(0.0100836),
    kd_(0.000206355),
    kT0_(273.15),
    kalpha_(5.0e-10),
    kp0_(1.0e5),
    M_(other.M_) {}

// ---------------------------------------------------------------------------
// Virtual copy constructor.
// ---------------------------------------------------------------------------
Teuchos::RCP<FieldModel> EOSWater::Clone() const {
  return Teuchos::rcp(new EOSWater(*this));
}

double EOSWater::Density(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_*dT)*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_));
};

double EOSWater::DDensityDT(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = kb_ + (2.0*kc_ + 3.0*kd_*dT)*dT;
  return rho1bar * (1.0 + kalpha_*(p - kp0_));

};

double EOSWater::DDensityDp(double T, double p) {
  double dT = T - kT0_;
  double rho1bar = ka_ + (kb_ + (kc_ + kd_*dT)*dT)*dT;
  return rho1bar * kalpha_;
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
