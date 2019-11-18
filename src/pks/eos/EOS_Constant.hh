/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Simple EOS for constant density.
  It defaults to reasonable values for water.
*/

#ifndef AMANZI_EOS_CONSTANT_HH_
#define AMANZI_EOS_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_ConstantMolarMass.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class EOS_Constant : public EOS_ConstantMolarMass {
 public:
  explicit EOS_Constant(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p) { return rho_; }
  virtual double DMassDensityDT(double T, double p) { return 0.0; }
  virtual double DMassDensityDp(double T, double p) { return 0.0; }

  virtual double DMolarDensityDT(double T, double p) { return 0.0; }
  virtual double DMolarDensityDp(double T, double p) { return 0.0; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;

  static Utils::RegisteredFactory<EOS, EOS_Constant> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
