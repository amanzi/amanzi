/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Quadratic internal energy model -- function of Cv and temperature

See ATS process model documentation's permafrost physical properties
documentation for details.

u = u0 + a(T - T_ref) + b(T - T_ref)^2 

UNITS: M/{mol,kg}
------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGYRELATIONS_IEM_QUADRATIC_
#define AMANZI_ENERGYRELATIONS_IEM_QUADRATIC_

#include "Teuchos_ParameterList.hpp"

#include "iem.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class IEMQuadratic : public IEM {

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
  static Utils::RegisteredFactory<IEM,IEMQuadratic> factory_;

};

}
}

#endif
