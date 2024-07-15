/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Thermal conductivity of ideal gas based on the Sutherland law.

         mu0 [W/m/K]  T0 [K]   S0 [K]

  Air    2.41e-2      273      194
  Argon  1.63e-2      273      170
  C02    1.46e-2      273     1800
  CO     2.32e-2      273      180
  N2     2.42e-2      273      240
  O2     2.44e-2      273      139
  H2     1.68e-2      273      120
  Steam  1.81e-2     300     2200

  Air: Rathakrishnan 2013: valid from 0.01 to 100 atm, and between 0 and 3000K. 
       Frank White’s "Viscous Fluid Flow":  2% error between 170K and 1900K.
*/

#ifndef AMANZI_EOS_IDEAL_GAS_THERMAL_CONDUCTIVITY_HH_
#define AMANZI_EOS_IDEAL_GAS_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

class IdealGas_ThermalConductivity : public EOS_ThermalConductivity {
 public:
  explicit IdealGas_ThermalConductivity(Teuchos::ParameterList& eos_plist);

  virtual double ThermalConductivity(double T, double p);
  virtual double DThermalConductivityDT(double T, double p);
  virtual double DThermalConductivityDP(double T, double p) { return 0.0; }

 private:
  double T0_, k0_, S0_, factor_;

 private:
  static Utils::RegisteredFactory<EOS_ThermalConductivity, IdealGas_ThermalConductivity> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
