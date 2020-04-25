/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for EOS implementations.
*/

#ifndef AMANZI_EOS_VAPOR_PRESSURE_BASE_FACTORY_HH_
#define AMANZI_EOS_VAPOR_PRESSURE_BASE_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VaporPressure_Base.hh"
#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

class VaporPressureBaseFactory : public Utils::Factory<VaporPressure_Base> {
 public:
  Teuchos::RCP<VaporPressure_Base> CreateVaporPressure(Teuchos::ParameterList& plist);
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
