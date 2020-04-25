/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Self-registering factory for flow PKs.
*/

#include "Darcy_PK.hh"
#include "Richards_PK.hh"
#include "VWContentEvaluator.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Darcy_PK> Darcy_PK::reg_("darcy");
RegisteredPKFactory<Richards_PK> Richards_PK::reg_("richards");
Utils::RegisteredFactory<FieldEvaluator, VWContentEvaluator> VWContentEvaluator::reg_("water content");

}  // namespace Amanzi
}  // namespace Flow

