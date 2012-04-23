/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Quadratic internal energy model -- function of Cv and temperature

See ATS process model documentation's permafrost physical properties
documentation for details.

u = u0 + a(T - T_ref) + b(T - T_ref)^2 

UNITS: J/{mol,kg}
------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGYRELATIONS_IEM_QUADRATIC_
#define AMANZI_ENERGYRELATIONS_IEM_QUADRATIC_

#include "Teuchos_ParameterList.hpp"

#include "internal_energy_model.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

class IEMQuadratic : public InternalEnergyModel {

public:
  explicit IEMQuadratic(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp);
  double DInternalEnergyDT(double temp);

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double u0_;
  double ka_;
  double kb_;
  double T0_; // units: K
  bool molar_basis_;

private:  
  // iem factor registration
  static Utils::RegisteredFactory<InternalEnergyModel,IEMQuadratic> factory_;

};

}
}
}

#endif
