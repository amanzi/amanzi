/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas (does not implement viscosity at this point!)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_IDEAL_GAS_HH_
#define _FLOWRELATIONS_EOS_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"
#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSIdealGas : public EOS {

public:
  explicit EOSIdealGas(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p);
  virtual double DMassDensityDT(double T, double p);
  virtual double DMassDensityDp(double T, double p);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

  virtual double molar_mass() { return M_; }

  virtual double Viscosity(double T) {}
  virtual double DViscosityDT(double T) {}

protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double R_;
  double M_;
};

} // namespace
} // namespace
} // namespace

#endif
