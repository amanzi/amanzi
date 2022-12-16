/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  Self-registering factory for thermal conductivity models.
*/

#ifndef PK_ENERGY_TCM_FACTORY_TWOPHASE_HH_
#define PK_ENERGY_TCM_FACTORY_TWOPHASE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TCM_TwoPhase.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class TCMFactory_TwoPhase : public Utils::Factory<TCM_TwoPhase> {
 public:
  Teuchos::RCP<TCM_TwoPhase> CreateTCM(Teuchos::ParameterList& plist);
};

} // namespace Energy
} // namespace Amanzi

#endif
