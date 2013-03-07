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
  return T < T0_ ? gamma_ * dens * (T0_ - T)/T0_ : 0.;
};

double PCIceWater::DCapillaryPressureDT(double T, double dens) {
  return T < T0_ ? -gamma_ * dens / T0_ : 0.;
};

double PCIceWater::DCapillaryPressureDRho(double T, double dens) {
  return T < T0_ ? gamma_ * (T0_ - T)/T0_ : 0.;
};


void PCIceWater::InitializeFromPlist_() {
  sigma_ice_liq_ = pc_plist_.get<double>("interfacial tension ice-water", 33.1);
  sigma_gas_liq_ = pc_plist_.get<double>("interfacial tension air-water", 72.7);
  T0_ = pc_plist_.get<double>("heat of fusion reference temperature [K]", 273.15);

  if (pc_plist_.isParameter("heat of fusion of water [J/mol]")) {
    heat_fusion_ = pc_plist_.get<double>("heat of fusion of water [J/mol]");
    molar_basis_ = true;
  } else {
    heat_fusion_ = pc_plist_.get<double>("heat of fusion of water [J/kg]", 3.34e5);
    molar_basis_ = false;
  }

  gamma_ = sigma_gas_liq_/sigma_ice_liq_ * heat_fusion_;
};

} // namespace
} // namespace
} // namespace
