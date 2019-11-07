/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for a combination of air and vapor pressure. Mass density is not
  available, not because it can't be calculated, but because it depends 
  upon omega. It's not really needed, and if it were, would not fit the 
  EOS interface without serious work.
*/

#ifndef AMANZI_EOS_VAPOR_IN_GAS_HH_
#define AMANZI_EOS_VAPOR_IN_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"

#include "EOS.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class EOS_VaporInGas : public EOS {
 public:
  EOS_VaporInGas(Teuchos::ParameterList& eos_plist);

  double MassDensity(double T, double p) { AMANZI_ASSERT(0); return 0.0; }
  double DMassDensityDT(double T, double p)  { AMANZI_ASSERT(0); return 0.0; }
  double DMassDensityDp(double T, double p)  { AMANZI_ASSERT(0); return 0.0; }

  double MolarDensity(double T, double p);
  double DMolarDensityDT(double T, double p);
  double DMolarDensityDp(double T, double p);

  bool IsConstantMolarMass() { return false; }
  double MolarMass() { AMANZI_ASSERT(0); return 0.0; }

protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  Teuchos::RCP<EOS> gas_eos_;

 private:
  static Utils::RegisteredFactory<EOS, EOS_VaporInGas> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
