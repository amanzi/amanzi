/*
  Multiphase PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Self-registering factory for multiphase PKs.
*/

#include "MultiphaseModelI_PK.hh"
#include "MultiphaseTwoComponents_PK.hh"

namespace Amanzi {
namespace Multiphase {

RegisteredPKFactory<MultiphaseModelI_PK> MultiphaseModelI_PK::reg_("multiphase pl sl xg");
RegisteredPKFactory<MultiphaseTwoComponents_PK> MultiphaseTwoComponents_PK::reg_("multiphase 2p2c");

}  // namespace Multiphase
}  // namespace Amanzi
