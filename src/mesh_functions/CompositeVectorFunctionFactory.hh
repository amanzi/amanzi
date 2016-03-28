/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Factory for a CV function on a mesh.
*/

#ifndef AMANZI_COMPOSITE_VECTOR_FUNCTION_FACTORY_HH_
#define AMANZI_COMPOSITE_VECTOR_FUNCTION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {
namespace Functions {

Teuchos::RCP<CompositeVectorFunction>
CreateCompositeVectorFunction(Teuchos::ParameterList& plist,
        const CompositeVectorSpace& factory);

} // namespace
} // namespace


#endif
