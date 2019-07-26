/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#include "Transport_PK.hh"
#include "TransportImplicit_PK.hh"

namespace Amanzi {
namespace Transport {

RegisteredPKFactory<Transport_PK> Transport_PK::reg_("transport");
RegisteredPKFactory<TransportImplicit_PK> TransportImplicit_PK::reg_("transport implicit");  

}  // namespace Transport
}  // namespace Amanzi
