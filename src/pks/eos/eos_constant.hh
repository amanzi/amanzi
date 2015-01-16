/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Simple EOS for constant density.
  Defaults to reasonable values for water.
*/

#ifndef AMANZI_EOS_CONSTANT_HH_
#define AMANZI_EOS_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos_constant_molar_mass.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class EOSConstant : public EOSConstantMolarMass {
 public:
  explicit EOSConstant(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p) { return rho_; }
  virtual double DMassDensityDT(double T, double p) { return 0.0; }
  virtual double DMassDensityDp(double T, double p) { return 0.0; }

  virtual double DMolarDensityDT(double T, double p) { return 0.0; }
  virtual double DMolarDensityDp(double T, double p) { return 0.0; }

 private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double rho_;

  static Utils::RegisteredFactory<EOS,EOSConstant> factory_;
};

}  // namespace Relations
}  // namespace Amanzi

#endif
