/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Model for effective permeability in fracutres.
*/

#ifndef AMANZI_FLOW_FRACTURE_PERM_MODEL_HH_
#define AMANZI_FLOW_FRACTURE_PERM_MODEL_HH_

#include <string>

namespace Amanzi {
namespace Flow {

class FracturePermModel {
 public:
  virtual ~FracturePermModel(){};
  virtual double Permeability(double aperture) = 0;
};

} // namespace Flow
} // namespace Amanzi

#endif
