/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "dbc.hh"
#include "pc_ice_water.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {


/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
PCIceWater::PCIceWater(Teuchos::ParameterList& pc_plist) :
    pc_plist_(pc_plist) {
  InitializeFromPlist_();
};

double PCIceWater::CapillaryPressure(double T, double dens) {
  return gamma_ * dens * (T0_ - T)/T0_;
};

double PCIceWater::DCapillaryPressureDT(double T, double dens) {
  return -gamma_ * dens / T0_;
};

double PCIceWater::DCapillaryPressureDRho(double T, double dens) {
  return -gamma_ * (T0_ - T)/T0_;
};


void PCIceWater::InitializeFromPlist_() {
  sigma_ice_liq_ = pc_plist_.get<double>("Interfacial tension ice-water");
  sigma_gas_liq_ = pc_plist_.get<double>("Interfacial tension air-water");
  T0_ = pc_plist_.get<double>("Heat of fusion reference temperature [K]", 273.15);

  if (pc_plist_.isParameter("Heat of fusion of water [J/kg]")) {
    ASSERT(0);  // not supported anymore!
    heat_fusion_ = pc_plist_.get<double>("Heat of fusion of water [J/kg]");
    molar_basis_ = false;
  } else {
    heat_fusion_ = pc_plist_.get<double>("Heat of fusion of water [J/mol]");
    molar_basis_ = true;
  }

  gamma_ = sigma_gas_liq_/sigma_ice_liq_ * heat_fusion_;
};

} // namespace
} // namespace
} // namespace
