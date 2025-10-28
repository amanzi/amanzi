/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  EOS for liquid water from package CoolProp
*/

#ifndef AMANZI_EOS_H2O_DENSITY_COOLPROP_HH_
#define AMANZI_EOS_H2O_DENSITY_COOLPROP_HH_

#include "AbstractState.h"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of state for water
class H2O_DensityCoolProp : public EOS_Density {
 public:
  H2O_DensityCoolProp(Teuchos::ParameterList& eos_plist);
  ~H2O_DensityCoolProp() { delete state_; }

  virtual double Density(double T, double p) override;
  virtual double DDensityDT(double T, double p) override;
  virtual double DDensityDp(double T, double p) override;

  virtual double MolarDensity(double T, double p) override { return Density(T, p) / M_; }
  virtual double DMolarDensityDT(double T, double p) override { return DDensityDT(T, p) / M_; }
  virtual double DMolarDensityDp(double T, double p) override { return DDensityDp(T, p) / M_; }

  CoolProp::phases get_phase(double T, double p);

 private:
  CoolProp::AbstractState* state_;

  static Utils::RegisteredFactory<EOS_Density, H2O_DensityCoolProp> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
