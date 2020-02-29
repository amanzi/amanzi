/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Linear internal energy model -- function of Cv and temperature

See ATS process model documentation's permafrost physical properties
documentation for details.

UNITS: MJ/{mol/kg}
------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGY_RELATIONS_IEM_LINEAR_
#define AMANZI_ENERGY_RELATIONS_IEM_LINEAR_

#include "Teuchos_ParameterList.hpp"

#include "iem.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class IEMLinear : public IEM {

public:
  explicit IEMLinear(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp);
  double DInternalEnergyDT(double temp) { return Cv_; }

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double Cv_; // units: MJ/({mol/kg}-K)
  double T_ref_; // units: K
  bool molar_basis_;
  double L_;

private:
  // iem factor registration
  static Utils::RegisteredFactory<IEM,IEMLinear> factory_;

};

} // namespace
} // namespace

#endif
