/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Factory for a CV function on a mesh.
------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITE_VECTOR_FUNCTION_FACTORY_HH_
#define AMANZI_COMPOSITE_VECTOR_FUNCTION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"

#include "composite_vector_function.hh"
#include "composite_vector_factory.hh"

namespace Amanzi {
namespace Functions {

Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        const CompositeVector& sample);

Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        const CompositeVectorFactory& factory);

} // namespace
} // namespace


#endif
