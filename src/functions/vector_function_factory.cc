/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $AMANZI_DIR/COPYRIGHT
Author Ethan Coon

Factory for functions from R^d to R^n.

See ATS_DIR/src/factory/factory.hh for an in-depth discussion of this
approach.

USAGE:

------------------------------------------------------------------------- */

#include "vector_function_factory.hh"

namespace Amanzi {

VectorFunctionFactory::map_type* VectorFunctionFactory::map_;

} // namespace
