/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear internal energy model -- function of Cv and temperature

See ATS process model documentation's permafrost physical properties
documentation for details.

UNITS: J/{mol/kg}
------------------------------------------------------------------------- */

#ifndef INTERNAL_ENERGY_LINEAR_
#define INTERNAL_ENERGY_LINEAR_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class InternalEnergyLinear {

public:
  InternalEnergyLinear(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp);
  double DInternalEnergyDT(double temp, double mol_frac_gas) { return Cv_; }

protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_; // units: J/({mol/kg}-K)
  double T_ref_; // units: K
  bool molar_basis_;
};

}
}
}

#endif
