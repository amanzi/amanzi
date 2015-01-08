/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  EOS for an ideal gas.
*/

#ifndef AMANZI_EOS_IDEAL_GAS_HH_
#define AMANZI_EOS_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSIdealGas : public EOSConstantMolarMass {
 public:
  explicit EOSIdealGas(Teuchos::ParameterList& eos_plist);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double R_;

 private:
  static Utils::RegisteredFactory<EOS,EOSIdealGas> factory_;
};

}  // namespace Relations
}  // namespace Amanzi

#endif
