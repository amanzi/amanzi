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

#ifndef AMANZI_EOS_DENSITY_VAPOR_IN_GAS_HH_
#define AMANZI_EOS_DENSITY_VAPOR_IN_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"

#include "EOS_Density.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class EOS_DensityVaporInGas : public EOS_Density {
 public:
  EOS_DensityVaporInGas(Teuchos::ParameterList& eos_plist);

  double Density(double T, double p) { AMANZI_ASSERT(0); return 0.0; }
  double DDensityDT(double T, double p)  { AMANZI_ASSERT(0); return 0.0; }
  double DDensityDp(double T, double p)  { AMANZI_ASSERT(0); return 0.0; }

  double MolarDensity(double T, double p) { return gas_eos_->MolarDensity(T, p); }
  double DMolarDensityDT(double T, double p) { return gas_eos_->DMolarDensityDT(T, p); }
  double DMolarDensityDp(double T, double p) { return gas_eos_->DMolarDensityDp(T, p); }

 protected:
  virtual void InitializeFromPlist_();

 protected:
  Teuchos::RCP<EOS_Density> gas_eos_;

 private:
  static Utils::RegisteredFactory<EOS_Density, EOS_DensityVaporInGas> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
