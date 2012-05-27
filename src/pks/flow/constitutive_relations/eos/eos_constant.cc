/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  Constant density/viscosity EOS, defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_constant.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<EOS,EOSConstant> EOSConstant::factory_("constant");

EOSConstant::EOSConstant(Teuchos::ParameterList& eos_plist) :
    eos_plist_(eos_plist) {

  InitializeFromPlist_();
};

void EOSConstant::InitializeFromPlist_() {
  // defaults to water
  if (eos_plist_.isParameter("Molar mass [kg/mol]")) {
    M_ = eos_plist_.get<double>("Molar mass [kg/mol]");
  } else {
    M_ = eos_plist_.get<double>("Molar mass [g/mol]", 18.0153) * 1.e-3;
  }

  rho_ = eos_plist_.get<double>("Density [kg/m^3]", 1000.0);
  visc_ = eos_plist_.get<double>("Viscosity [Pa-s]", 8.9e-4);

};

} // namespace
} // namespace
} // namespace
