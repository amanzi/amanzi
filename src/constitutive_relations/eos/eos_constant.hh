/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  Simple EOS for constant density.
  Defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_CONSTANT_HH_
#define AMANZI_RELATIONS_EOS_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSConstant : public EOSConstantMolarMass {

public:
  explicit EOSConstant(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p) { return rho_; }
  virtual double DMassDensityDT(double T, double p) { return 0.0; }
  virtual double DMassDensityDp(double T, double p) { return 0.0; }

  virtual double DMolarDensityDT(double T, double p) { return 0.0; }
  virtual double DMolarDensityDp(double T, double p) { return 0.0; }

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;

  static Utils::RegisteredFactory<EOS,EOSConstant> factory_;

};

} // namespace
} // namespace

#endif
