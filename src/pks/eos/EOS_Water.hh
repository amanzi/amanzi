/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water for T between 0 and 30 C. Need a reference.
*/

#ifndef AMANZI_EOS_LIQUID_WATER_HH_
#define AMANZI_EOS_LIQUID_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EOS_ConstantMolarMass.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class EOS_Water : public EOS_ConstantMolarMass {
 public:
  explicit EOS_Water(Teuchos::ParameterList& eos_plist);

  virtual double MassDensity(double T, double p);
  virtual double DMassDensityDT(double T, double p);
  virtual double DMassDensityDp(double T, double p);

 private:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  const double ka_, kb_, kc_, kd_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  static Utils::RegisteredFactory<EOS, EOS_Water> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
