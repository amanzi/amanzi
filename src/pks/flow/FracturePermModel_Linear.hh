/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Linear model for effective permeability in fracutres.
*/

#ifndef AMANZI_FRACTURE_PERM_MODEL_LINEAR_HH_
#define AMANZI_FRACTURE_PERM_MODEL_LINEAR_HH_

#include "Teuchos_ParameterList.hpp"

#include "FracturePermModel.hh"

namespace Amanzi {
namespace Flow {

class FracturePermModel_Linear : public FracturePermModel {
 public:
  explicit FracturePermModel_Linear(Teuchos::ParameterList& plist) {}
  ~FracturePermModel_Linear(){};

  // required methods from the base class
  inline double Permeability(double aperture) { return aperture / 12; }
};

} // namespace Flow
} // namespace Amanzi

#endif
