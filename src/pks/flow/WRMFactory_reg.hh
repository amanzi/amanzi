/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for WRM implementations.
*/

#include "WRMFactory.hh"

// explicity instantitate the static data of factory
namespace Amanzi {
namespace Utils {

template<> Factory<Flow::WRM>::map_type* Factory<Flow::WRM>::map_;

}  // namespace Utils
}  // namespace Amanzi

