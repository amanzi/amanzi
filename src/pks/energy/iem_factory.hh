/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Self-registering factory for IEM implementations.
*/

#ifndef AMANZI_ENERGYRELATIONS_IEM_FACTORY_
#define AMANZI_ENERGYRELATIONS_IEM_FACTORY_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "iem.hh"
#include "factory.hh"

namespace Amanzi {
namespace Energy {

class IEMFactory : public Utils::Factory<IEM> {
 public:
  Teuchos::RCP<IEM> createIEM(Teuchos::ParameterList& plist);
};

}  // namespace Energy
}  // namespace Amanzi

#endif
