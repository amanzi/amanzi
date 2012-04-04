/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_WATER_HH_
#define _FLOWRELATIONS_EOS_WATER_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSWater {

public:
  explicit EOSWater(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p);
  virtual double DMassDensityDT(double T, double p);
  virtual double DMassDensityDp(double T, double p);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

  virtual double Viscosity(double T);

  double molar_mass() { return M_; }

protected:
  virtual void InitializeFromPlist_();
  virtual double DViscosityDT(double T); // not yet implemented

  Teuchos::ParameterList eos_plist_;
  double M_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  const double ka_, kb_, kc_, kd_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  // -- temperature dependence of viscosity < T1
  const double kav1_, kbv1_, kcv1_;

  // -- temperature dependence of viscosity > T1
  const double kbv2_, kcv2_, kT1_;
};

} // namespace
} // namespace
} // namespace

#endif
