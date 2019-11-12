/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Self-registering factory for MDM implementations.
*/

#ifndef PK_TRANSPORT_MDM_FACTORY_HH_
#define PK_TRANSPORT_MDM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MDM.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Transport {

class MDMFactory : public Utils::Factory<MDM> {
 public:
  Teuchos::RCP<MDM> Create(Teuchos::ParameterList& plist);
};

}  // namespace Transport
}  // namespace Amanzi

#endif
