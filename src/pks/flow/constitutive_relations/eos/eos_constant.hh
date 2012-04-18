/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  Simple EOS for constant density and viscosity.
  Defaults to reasonable values for water.

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_CONSTANT_HH_
#define _FLOWRELATIONS_EOS_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSConstant : public EOS {

public:
  explicit EOSConstant(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p) { return rho_; }
  virtual double DMassDensityDT(double T, double p) { return 0.0; }
  virtual double DMassDensityDp(double T, double p) { return 0.0; }

  virtual double MolarDensity(double T, double p) { return rho_/M_; }
  virtual double DMolarDensityDT(double T, double p) { return 0.0; }
  virtual double DMolarDensityDp(double T, double p) { return 0.0; }

  virtual double Viscosity(double T) { return visc_; }

  double molar_mass() { return M_; }

  double DViscosityDT(double T) { return 0.0; } // not implemented

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double M_;
  double rho_;
  double visc_;

  static Utils::RegisteredFactory<EOS,EOSConstant> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
