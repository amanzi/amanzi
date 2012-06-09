/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Basics of an EOS class... likely isn't specific enough for any given EOS,
  and will need other base models.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _PK_FLOW_EOS_HH_
#define _PK_FLOW_EOS_HH_

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOS {

public:
  // constructor format for all derived classes
  // explicit EOS(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p) = 0;
  virtual double DMassDensityDT(double T, double p) = 0;
  virtual double DMassDensityDp(double T, double p) = 0;

  virtual double MolarDensity(double T, double p) = 0;
  virtual double DMolarDensityDT(double T, double p) = 0;
  virtual double DMolarDensityDp(double T, double p) = 0;

  virtual double molar_mass() = 0;

  virtual double Viscosity(double T) = 0;
  virtual double DViscosityDT(double T) = 0;

};

} // namespace
} // namespace
} // namespace

#endif
