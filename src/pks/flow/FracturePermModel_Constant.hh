/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Constant model for effective permeability in a fracture.
*/

#ifndef AMANZI_FRACTURE_PERM_MODEL_CONSTANT_HH_
#define AMANZI_FRACTURE_PERM_MODEL_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "FracturePermModel.hh"

namespace Amanzi {
namespace Flow {

class FracturePermModel_Constant : public FracturePermModel {
 public:
  explicit FracturePermModel_Constant(Teuchos::ParameterList& plist)
  {
    value_ = plist.get<double>("value");
  }
  ~FracturePermModel_Constant(){};

  // required methods from the base class
  inline double Permeability(double aperture) { return value_; }

 private:
  double value_;
};

} // namespace Flow
} // namespace Amanzi

#endif
