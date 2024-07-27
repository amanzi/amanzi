/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

/*
 Shallow water PK

 */

#include "WhetStoneDefs.hh"

#include "MechanicsBoundaryFunction.hh"

namespace Amanzi {
namespace Mechanics {

/* ****************************************************************
* Constructor: extract attributes to setup a submodel.
**************************************************************** */
MechanicsBoundaryFunction::MechanicsBoundaryFunction(const Teuchos::ParameterList& plist)
{
  if (plist.isParameter("plane strain direction"))
    plane_strain_direction_ = plist.get<std::string>("plane strain direction");
}

} // namespace Mechanics
} // namespace Amanzi
