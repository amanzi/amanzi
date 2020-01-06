/*
  Multiphase PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Self-registering factory for multiphase PKs.
*/

#include "MultiphaseReduced_PK.hh"

namespace Amanzi {
namespace Multiphase {

RegisteredPKFactory<MultiphaseReduced_PK> MultiphaseReduced_PK::reg_("multiphase reduced");

}  // namespace Multiphase
}  // namespace Amanzi
