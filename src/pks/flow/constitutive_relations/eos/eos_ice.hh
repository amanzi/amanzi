/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid ice.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_ICE_HH_
#define _FLOWRELATIONS_EOS_ICE_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSIce : public EOS {

public:
  explicit EOSIce(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p);
  virtual double DMassDensityDT(double T, double p);
  virtual double DMassDensityDp(double T, double p);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

  double Viscosity(double T) {}

  double molar_mass() { return M_; }

  double DViscosityDT(double T) {} // not implemented

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double M_;

  // constants for ice, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  const double ka_, kb_, kc_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  static Utils::RegisteredFactory<EOS,EOSIce> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
