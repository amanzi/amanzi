/*
  Multiphase PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Self-registering factory for multiphase PKs.
*/

#include "MultiphaseJaffre_PK.hh"
#include "MultiphaseModel1_PK.hh"
#include "MultiphaseModel2_PK.hh"

namespace Amanzi {
namespace Multiphase {

RegisteredPKFactory<MultiphaseJaffre_PK> MultiphaseJaffre_PK::reg_("multiphase jaffre");
RegisteredPKFactory<MultiphaseModel1_PK> MultiphaseModel1_PK::reg_("multiphase pl sl xg");
RegisteredPKFactory<MultiphaseModel2_PK> MultiphaseModel2_PK::reg_("multiphase pl sl ng");

}  // namespace Multiphase
}  // namespace Amanzi
