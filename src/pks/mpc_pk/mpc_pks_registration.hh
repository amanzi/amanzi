/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

 FLowReactiveTransport_PK registration
*/

#include "FlowEnergy_PK.hh"
#include "FlowMatrixFracture_PK.hh"
#include "FlowReactiveTransport_PK.hh"
#include "ReactiveTransport_PK.hh"

namespace Amanzi {

RegisteredPKFactory<FlowMatrixFracture_PK> FlowMatrixFracture_PK::reg_("darcy matrix fracture");
RegisteredPKFactory<FlowEnergy_PK> FlowEnergy_PK::reg_("thermal richards");
RegisteredPKFactory<FlowReactiveTransport_PK> FlowReactiveTransport_PK::reg_("flow reactive transport");
RegisteredPKFactory<ReactiveTransport_PK> ReactiveTransport_PK::reg_("reactive transport");

}  // namespace Amanzi

