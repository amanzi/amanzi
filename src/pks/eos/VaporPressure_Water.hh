/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Saturated vapor pressure for vapor over water or ice, Sonntag (1990)
*/

#ifndef AMANZI_EOS_VAPOR_PRESSURE_WATER_HH_
#define AMANZI_EOS_VAPOR_PRESSURE_WATER_HH_

#include "Teuchos_ParameterList.hpp"
#include "Factory.hh"
#include "VaporPressure_Base.hh"

namespace Amanzi {
namespace AmanziEOS {

class VaporPressure_Water : public VaporPressure_Base {
 public:
  explicit VaporPressure_Water(Teuchos::ParameterList& plist);

  virtual double SaturatedVaporPressure(double T);
  virtual double DSaturatedVaporPressureDT(double T);

 private:
  Teuchos::ParameterList plist_;
  const double ka0_;
  const double ka_, kb_, kc_, kd_;

  static Utils::RegisteredFactory<VaporPressure_Base, VaporPressure_Water> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
