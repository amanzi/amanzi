/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Cubic law model for effective permeability in fracutres.
*/

#ifndef AMANZI_FRACTURE_PERM_MODEL_CUBICLAW_HH_
#define AMANZI_FRACTURE_PERM_MODEL_CUBICLAW_HH_

#include "Teuchos_ParameterList.hpp"

#include "FracturePermModel.hh"

namespace Amanzi {
namespace Flow {

class FracturePermModel_CubicLaw : public FracturePermModel {
 public:
  explicit FracturePermModel_CubicLaw(Teuchos::ParameterList& plist){};
  ~FracturePermModel_CubicLaw(){};

  // required methods from the base class
  inline double Permeability(double aperture) { return aperture * aperture / 12; }
};

} // namespace Flow
} // namespace Amanzi

#endif
