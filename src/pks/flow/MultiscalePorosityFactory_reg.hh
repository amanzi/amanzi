/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Self-registering factory for multiscale porosity models.
*/

#include "MultiscalePorosityFactory.hh"

// explicity instantitate the static data of factory
namespace Amanzi {
namespace Utils {

template<>
Factory<Flow::MultiscalePorosity>::map_type* Factory<Flow::MultiscalePorosity>::map_;

}  // namespace Utils
}  // namespace Amanzi

