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

#ifndef AMANZI_ENERGYRELATIONS_IEM_LINEAR_
#define AMANZI_ENERGYRELATIONS_IEM_LINEAR_

#include "Teuchos_ParameterList.hpp"

#include "internal_energy_model.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class IEMLinear : public InternalEnergyModel {

public:
  explicit IEMLinear(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp);
  double DInternalEnergyDT(double temp) { return Cv_; }

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_; // units: J/({mol/kg}-K)
  double T_ref_; // units: K
  bool molar_basis_;

private:  
  // iem factor registration
  static Utils::RegisteredFactory<InternalEnergyModel,IEMLinear> factory_;

};

}
}
}

#endif
