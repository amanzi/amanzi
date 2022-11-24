/*
   Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
   Amanzi is released under the three-clause BSD License.
   The terms of use and "as is" disclaimer for this license are
   provided in the top-level COPYRIGHT file.
   See $ATS_DIR/COPYRIGHT

  Author: Ethan Coon (ecoon@lanl.gov)
*/

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "MeshFunction.hh"
#include "CompositeVector.hh"

namespace Amanzi {

class VerboseObject;

namespace Functions {

class CompositeVectorFunction final {
 public:
  CompositeVectorFunction(const Teuchos::RCP<const MeshFunction>& func,
                          const std::vector<std::string>& names);
  ~CompositeVectorFunction() = default;

  void
  Compute(double time, const Teuchos::Ptr<CompositeVector>& vec, const VerboseObject* vo = nullptr);

 protected:
  typedef std::pair<std::string, Teuchos::RCP<MeshFunction::Spec>> CompositeVectorSpec;
  typedef std::vector<Teuchos::RCP<CompositeVectorSpec>> CompositeVectorSpecList;

  Teuchos::RCP<const MeshFunction> func_;
  CompositeVectorSpecList cv_spec_list_;
};

} // namespace Functions
} // namespace Amanzi
