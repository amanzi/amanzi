/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Simple EOS for compressibility with pressure.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_EOS_LINEAR_HH_
#define AMANZI_RELATIONS_EOS_LINEAR_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSLinear : public EOSConstantMolarMass {

public:
  explicit EOSLinear(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p) { return rho_ * (1+beta_*std::max(p-101325., 0.)); }
  virtual double DMassDensityDp(double T, double p) { return p > 101325. ? rho_ * beta_ : 0.; }
  virtual double DMassDensityDT(double T, double p) { return 0.; }

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;
  double beta_;

  static Utils::RegisteredFactory<EOS,EOSLinear> factory_;

};

} // namespace
} // namespace

#endif
