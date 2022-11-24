/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_POROSITY_MODEL_CONSTANT_HH_
#define AMANZI_POROSITY_MODEL_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "PorosityModel.hh"

namespace Amanzi {
namespace Flow {

class PorosityModel_Constant : public PorosityModel {
 public:
  explicit PorosityModel_Constant(Teuchos::ParameterList& plist)
  {
    porosity_ = plist.get<double>("value");
  }
  ~PorosityModel_Constant(){};

  // required methods from the base class
  inline double Porosity(double p) { return porosity_; }
  inline double dPorositydPressure(double p) { return 0.0; }

 private:
  double porosity_;
};

} // namespace Flow
} // namespace Amanzi

#endif
