/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#ifndef AMANZI_COMPOSITE_VECTOR_FUNCTION_HH_
#define AMANZI_COMPOSITE_VECTOR_FUNCTION_HH_

#include <string>
#include <utility>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "MeshFunction.hh"
#include "CompositeVector.hh"

namespace Amanzi {
namespace Functions {

class CompositeVectorFunction {
 public:
  CompositeVectorFunction(const Teuchos::RCP<const MeshFunction>& func,
                          const std::vector<std::string>& names);
  virtual ~CompositeVectorFunction() = default;

  virtual void Compute(double time, CompositeVector& vec);

 protected:
  typedef std::pair<std::string, Teuchos::RCP<MeshFunction::Spec>>
    CompositeVectorSpec;
  typedef std::vector<Teuchos::RCP<CompositeVectorSpec>>
    CompositeVectorSpecList;

  Teuchos::RCP<const MeshFunction> func_;
  CompositeVectorSpecList cv_spec_list_;
};

} // namespace Functions
} // namespace Amanzi

#endif
