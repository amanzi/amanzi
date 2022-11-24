/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Saturated vapor pressure for vapor over water or ice, Sonntag (1990)
*/

#ifndef AMANZI_EOS_H2O_SATURATED_VAPOR_PRESSURE_HH_
#define AMANZI_EOS_H2O_SATURATED_VAPOR_PRESSURE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "EOS_SaturatedVaporPressure.hh"

namespace Amanzi {
namespace AmanziEOS {

class H2O_SaturatedVaporPressure : public EOS_SaturatedVaporPressure {
 public:
  explicit H2O_SaturatedVaporPressure(Teuchos::ParameterList& plist);

  virtual double Pressure(double T);
  virtual double DPressureDT(double T);

 private:
  Teuchos::ParameterList plist_;
  const double ka0_;
  const double ka_, kb_, kc_, kd_;

  static Utils::RegisteredFactory<EOS_SaturatedVaporPressure, H2O_SaturatedVaporPressure> factory_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
