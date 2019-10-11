/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/
//! Mesh Functions, evaluate a function on a mesh and stick the result in a
//! vector.

/*!

CompositeVectorFunctions are ways of evaluating a piecewise function on a
mesh and sticking the result into a CompositeVector.

This is used in a variety of ways -- Initial Conditions, Boundary
Conditions, and Independent Variable Evaluators.

Typically any containing object of this spec is a list of these specs.  The
list is indexed by Region, and the regions (logically) should partition the
domain (or boundary of the domain in the case of BCs).

Each entry in that list is a:

``[composite-vector-function-spec]``

ONE OF:
* `"region`" ``[string]`` Region on which this function is evaluated.
OR:
* `"regions`" ``[Array(string)]`` List of regions on which this function is
evaluated. END

ONE OF:
* `"component`" ``[string]`` Mesh component to evaluate this on.  This is one of
"cell", "face", or "node" OR:
* `"components`" ``[Array(string)]`` Mesh components to evaluate this on.  This
is some collection of "cell", "face", and/or "node" END

* `"function`" ``[function-spec]`` The spec to provide the actual algebraic
function.

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

} // namespace Functions
} // namespace Amanzi


#endif
