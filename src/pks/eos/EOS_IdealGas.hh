/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for the ideal gas.
*/

#ifndef AMANZI_EOS_IDEAL_GAS_HH_
#define AMANZI_EOS_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_ConstantMolarMass.hh"

namespace Amanzi {
namespace AmanziEOS {
 
// Equation of State model
class EOS_IdealGas : public EOS_ConstantMolarMass {
 public:
  explicit EOS_IdealGas(Teuchos::ParameterList& eos_plist);

  virtual double MolarDensity(double T, double p);
  virtual double DMolarDensityDT(double T, double p);
  virtual double DMolarDensityDp(double T, double p);

 protected:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double R_;

 private:
  static Utils::RegisteredFactory<EOS, EOS_IdealGas> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
