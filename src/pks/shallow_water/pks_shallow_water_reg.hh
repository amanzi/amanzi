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


#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

RegisteredPKFactory<ShallowWater_PK> ShallowWater_PK::reg_("shallow water");

} // namespace ShallowWater
} // namespace Amanzi
