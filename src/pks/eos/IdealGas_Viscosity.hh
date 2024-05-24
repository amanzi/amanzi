/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Viscosity for ideal gas based on the Sutherland law.

         mu0 [Pa.s]   T0 [K]   S0 [K]

  Air    1.716e-5     273      111
  Argon  2.125e-5     273      114
  C02    1.370e-5     273      222
  CO     1.657e-5     273      136
  N2     1.663e-5     273      107
  O2     1.919e-5     273      139
  H2     8.411e-5     273       97
  Steam  1.120e-5     350     1064

  Air: Rathakrishnan 2013: valid from 0.01 to 100 atm, and between 0 and 3000K. 
       Frank White’s "Viscous Fluid Flow":  2% error between 170K and 1900K.
*/

#ifndef AMANZI_EOS_IDEAL_GAS_VISCOSITY_HH_
#define AMANZI_EOS_IDEAL_GAS_VISCOSITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

class IdealGas_Viscosity : public EOS_Viscosity {
 public:
  explicit IdealGas_Viscosity(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T, double p);
  virtual double DViscosityDT(double T, double p);
  virtual double DViscosityDp(double T, double p) { return 0.0; }

 private:
  double T0_, mu0_, S0_, factor_;

 private:
  static Utils::RegisteredFactory<EOS_Viscosity, IdealGas_Viscosity> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
