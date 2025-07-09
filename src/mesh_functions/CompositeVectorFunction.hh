/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
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
                          const std::vector<std::string>& names,
                          bool dot_with_normal = false,
                          const std::string& spatial_dist_method = "none");
  ~CompositeVectorFunction() = default;

  void
  Compute(double time, const Teuchos::Ptr<CompositeVector>& vec, const VerboseObject* vo = nullptr);

 protected:
  void
  Compute_(double time, const Teuchos::Ptr<CompositeVector>& vec, const VerboseObject* vo);

  void
  ComputeDotWithNormal_(double time, const Teuchos::Ptr<CompositeVector>& vec, const VerboseObject* vo);

  void
  ComputeSpatiallyDistributed_(double time, const Teuchos::Ptr<CompositeVector>& vec, const VerboseObject* vo);

  // compute a single spec
  void
  ComputeSpec_(const MeshFunction::Spec& spec,
               double time,
               Epetra_MultiVector& vec,
               const VerboseObject* vo);

 protected:
  typedef std::pair<std::string, Teuchos::RCP<MeshFunction::Spec>> CompositeVectorSpec;
  typedef std::vector<Teuchos::RCP<CompositeVectorSpec>> CompositeVectorSpecList;

  Teuchos::RCP<const MeshFunction> func_;
  CompositeVectorSpecList cv_spec_list_;
  std::string spatial_dist_method_;
  bool dot_with_normal_;
};

} // namespace Functions
} // namespace Amanzi
