/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include "FlowEnergy_PK.hh"

namespace Amanzi {

RegisteredPKFactory<FlowEnergy_PK> FlowEnergy_PK::reg_("flow energy");

}  // namespace Amanzi
