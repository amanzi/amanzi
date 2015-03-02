/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for WRM implementations.
*/

#include "WRMFactory.hh"

namespace Amanzi {
namespace Flow {

// explicity instantitate the static data of factory
template<> 
Amanzi::Utils::Factory<Amanzi::Flow::WRM>::map_type* 
Amanzi::Utils::Factory<Amanzi::Flow::WRM>::map_;

}  // namespace Flow
}  // namespace Amanzi

