/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Linear internal energy model is function of cv and temperature.
  UNITS: J/{mol/kg}
*/

#ifndef AMANZI_ENERGY_IEM_LINEAR_HH_
#define AMANZI_ENERGY_IEM_LINEAR_HH_

#include "Teuchos_ParameterList.hpp"

#include "IEM.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {

class IEM_Linear : public IEM {
 public:
  explicit IEM_Linear(Teuchos::ParameterList& plist);

  bool IsMolarBasis() { return molar_basis_; }

  double InternalEnergy(double temp);
  double DInternalEnergyDT(double temp) { return cv_; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double cv_;  // units: J/({mol/kg}-K)
  double Tref_;  // units: K
  bool molar_basis_;

 private:
  // iem factor registration
  static Utils::RegisteredFactory<IEM,IEM_Linear> factory_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
