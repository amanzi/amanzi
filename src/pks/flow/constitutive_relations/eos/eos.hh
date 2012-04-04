/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Basics of an EOS class... likely isn't specific enough for any given EOS,
  and will need other base models.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _EOS_HH_
#define _EOS_HH_

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOS {

public:
  explicit EOS(Teuchos::ParameterList& eos_plist);
  virtual double MassDensity(double T, double p) = 0;
  virtual double MolarDensity(double T, double p) = 0;
  virtual double Viscosity(double T) = 0;
private:
  void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
};

} // namespace
} // namespace
} // namespace

#endif
