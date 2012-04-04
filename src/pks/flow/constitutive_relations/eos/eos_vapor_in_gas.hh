/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas air with a molar fraction of water vapor.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_VAPOR_IN_GAS_HH_
#define _FLOWRELATIONS_EOS_VAPOR_IN_GAS_HH_

#include "Teuchos_ParameterList.hpp"
#include "eos_ideal_gas.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSVaporInGas : public EOSIdealGas {

public:
  explicit EOSVaporInGas(Teuchos::ParameterList& eos_plist) : EOSIdealGas(eos_plist) {}

  virtual double MassDensity(double T, double p, double mol_frac_vapor);
  virtual double DMassDensityDT(double T, double p, double mol_frac_vapor);
  virtual double DMassDensityDp(double T, double p, double mol_frac_vapor);
  virtual double DMassDensityDmol_frac(double T, double p, double mol_frac_vapor);

  virtual double MolarMass(double mol_frac_vapor) {
    return mol_frac_vapor*Mv_ + (1.0-mol_frac_vapor)*Mg_; }

  virtual double SaturatedVaporPressure(double T);

protected:
  virtual double MolarMass(); // undefined
  virtual double MassDensity(double T, double p);  // undefined, requires mol fraction
  virtual double Viscosity(double T); // undefined at this point -- make public when defined

  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double Mg_;
  double Mv_;
};

} // namespace
} // namespace
} // namespace

#endif
