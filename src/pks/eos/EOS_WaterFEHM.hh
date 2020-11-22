/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water for T between 0.001 and 360 C from FEHM manual
*/

#ifndef AMANZI_EOS_LIQUID_WATER_FEHM_HH_
#define AMANZI_EOS_LIQUID_WATER_FEHM_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_ConstantMolarMass.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class EOS_WaterFEHM : public EOS_ConstantMolarMass {
 public:
  explicit EOS_WaterFEHM(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p);
  virtual double DMassDensityDT(double T, double p);
  virtual double DMassDensityDp(double T, double p);

 private:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  double y0_, y1_, y2_, y3_, y4_, y5_, y6_, y7_, y8_, y9_;
  double z0_, z1_, z2_, z3_, z4_, z5_, z6_, z7_, z8_, z9_;
  double T0_;

  static Utils::RegisteredFactory<EOS, EOS_WaterFEHM> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
