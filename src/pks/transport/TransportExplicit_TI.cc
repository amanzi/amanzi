/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include <algorithm>

#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "TransportExplicit_PK.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Routine takes single component and returns functional value
****************************************************************** */
void
TransportExplicit_PK::FunctionalTimeDerivative(double t,
                                               const CompositeVector& component,
                                               CompositeVector& f)
{
  if (method_ == Method_t::MUSCL)
    FunctionalTimeDerivative_MUSCL_(t, component, f, true);
  else
    FunctionalTimeDerivative_FCT_(t, component, f);
}

} // namespace Transport
} // namespace Amanzi
