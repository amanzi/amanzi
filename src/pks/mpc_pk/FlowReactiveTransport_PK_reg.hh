/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

 FLowReactiveTransport_PK registration
*/

#include "FlowReactiveTransport_PK.hh"

namespace Amanzi {

RegisteredPKFactory<FlowReactiveTransport_PK> FlowReactiveTransport_PK::reg_("flow reactive transport");

}  // namespace Amanzi
