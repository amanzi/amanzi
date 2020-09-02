/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Saturated vapor pressure for vapor over water or ice, Sonntag (1990)
*/

#ifndef AMANZI_EOS_SATURATED_VAPOR_PRESSURE_WATER_HH_
#define AMANZI_EOS_SATURATED_VAPOR_PRESSURE_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "SaturatedVaporPressure.hh"

namespace Amanzi {
namespace AmanziEOS {

class SaturatedVaporPressure_Water : public SaturatedVaporPressure {
 public:
  explicit SaturatedVaporPressure_Water(Teuchos::ParameterList& plist);

  virtual double Pressure(double T);
  virtual double DPressureDT(double T);

 private:
  Teuchos::ParameterList plist_;
  const double ka0_;
  const double ka_, kb_, kc_, kd_;

  static Utils::RegisteredFactory<SaturatedVaporPressure, SaturatedVaporPressure_Water> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
