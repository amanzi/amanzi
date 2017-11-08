/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  A capillary pressure model based upon something other than p_atm - p.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_PC_ICE_WATER_
#define AMANZI_FLOW_RELATIONS_PC_ICE_WATER_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class PCIceWater {

public:
  explicit PCIceWater(Teuchos::ParameterList& plist);

  // required methods from the base class
  bool IsMolarBasis() { return molar_basis_; }
  double CapillaryPressure(double T, double dens);
  double DCapillaryPressureDT(double T, double dens);
  double DCapillaryPressureDRho(double T, double dens);

private:
  void InitializeFromPlist_();

  Teuchos::ParameterList pc_plist_;

  bool molar_basis_;
  double sigma_ice_liq_;
  double sigma_gas_liq_;
  double heat_fusion_;
  double gamma_;
  double T0_;
  double halfwidth_;
};

} //namespace
} //namespace

#endif
