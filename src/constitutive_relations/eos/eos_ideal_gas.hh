/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_IDEAL_GAS_HH_
#define AMANZI_RELATIONS_EOS_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSIdealGas : public EOSConstantMolarMass {

public:
  explicit EOSIdealGas(Teuchos::ParameterList& eos_plist);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double R_;

private:
  static Utils::RegisteredFactory<EOS,EOSIdealGas> factory_;
};

} // namespace
} // namespace

#endif
