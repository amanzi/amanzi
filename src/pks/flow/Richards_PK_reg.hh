/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Richards pprocess kernel registration.
*/

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards_PK> Richards_PK::reg_("richards");

}  // namespace Flow
}  // namespace Amanzi
