/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  Self-registering factory for IEM implementations.
*/

#ifndef AMANZI_ENERGY_IEM_FACTORY_HH_
#define AMANZI_ENERGY_IEM_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "IEM.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class IEMFactory : public Utils::Factory<IEM> {
 public:
  Teuchos::RCP<IEM> CreateIEM(Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace Amanzi

#endif
