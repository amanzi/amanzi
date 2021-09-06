/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Base class for saturated vapor pressure.
*/

#ifndef AMANZI_EOS_SATURATED_VAPOR_PRESSURE_HH_
#define AMANZI_EOS_SATURATED_VAPOR_PRESSURE_HH_

namespace Amanzi {
namespace AmanziEOS {

class EOS_SaturatedVaporPressure {
 public:
  virtual ~EOS_SaturatedVaporPressure() {};

  virtual double Pressure(double T) = 0;
  virtual double DPressureDT(double T) = 0;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
