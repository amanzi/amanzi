/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Function applied to a mesh component, along with meta-data to store the values
of this function in a ComposteVector.

------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITE_VECTOR_FUNCTION_HH_
#define AMANZI_COMPOSITE_VECTOR_FUNCTION_HH_

#include <string>
#include <utility>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "mesh_function.hh"
#include "composite_vector.hh"

namespace Amanzi {
namespace Functions {

class CompositeVectorFunction {

public:
  CompositeVectorFunction(const Teuchos::RCP<const MeshFunction>& func,
                          const std::vector<std::string>& names);

  virtual void Compute(double time, const Teuchos::Ptr<CompositeVector>& vec);

protected:
  typedef std::pair<std::string, Teuchos::RCP<MeshFunction::Spec> > CompositeVectorSpec;
  typedef std::vector<Teuchos::RCP<CompositeVectorSpec> > CompositeVectorSpecList;

  Teuchos::RCP<const MeshFunction> func_;
  CompositeVectorSpecList cv_spec_list_;
};

} // namespace
} // namespace

#endif
