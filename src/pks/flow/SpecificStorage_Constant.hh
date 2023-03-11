/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_FLOW_SPECIFIC_STORAGE_CONSTANT_MODEL_HH_
#define AMANZI_FLOW_SPECIFIC_STORAGE_CONSTANT_MODEL_HH_

#include <string>

#include "Factory.hh"

#include "SpecificStorage.hh"

namespace Amanzi {
namespace Flow {

class SpecificStorage_Constant : public SpecificStorage {
 public:
  SpecificStorage_Constant(Teuchos::ParameterList& plist) : SpecificStorage(plist)
  {
    value_ = plist.get<double>("value");
  }

  virtual double Value(double porosity) const { return value_; }

 private:
  double value_;

  static Utils::RegisteredFactory<SpecificStorage, SpecificStorage_Constant> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
