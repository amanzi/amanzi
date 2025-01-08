/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

/*
  Pipe flow PK

*/


#include "PipeFlow_PK.hh"

namespace Amanzi {
namespace ShallowWater {

RegisteredPKFactory<PipeFlow_PK> PipeFlow_PK::reg_("pipe flow");

} // namespace ShallowWater
} // namespace Amanzi
