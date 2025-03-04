/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_FLOW_SPECIFIC_STORAGE_MODEL_HH_
#define AMANZI_FLOW_SPECIFIC_STORAGE_MODEL_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class SpecificStorage {
 public:
  SpecificStorage(Teuchos::ParameterList& plist){};
  virtual ~SpecificStorage(){};

  virtual double Value(double porosity) const = 0;
};

} // namespace Flow
} // namespace Amanzi

#endif
