/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_POROSITY_MODEL_HH_
#define AMANZI_POROSITY_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class PorosityModel {
 public:
  virtual ~PorosityModel(){};
  virtual double Porosity(double p) = 0;
  virtual double dPorositydPressure(double p) = 0; // derivative wrt to pressure
};

} // namespace Flow
} // namespace Amanzi

#endif
