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

#ifndef AMANZI_EOS_H2O_VISCOSITY_COOLPROP_HH_
#define AMANZI_EOS_H2O_VISCOSITY_COOLPROP_HH_

#include "AbstractState.h"
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of state for water
class H2O_ViscosityCoolProp : public EOS_Viscosity {
 public:
  H2O_ViscosityCoolProp(Teuchos::ParameterList& eos_plist);
  ~H2O_ViscosityCoolProp() { delete state_; }

  virtual double Viscosity(double T, double p) override;
  virtual double DViscosityDT(double T, double p) override;
  virtual double DViscosityDp(double T, double p) override;

  CoolProp::phases get_phase(double T, double p);

 private:
  CoolProp::AbstractState* state_;

  static Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityCoolProp> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
