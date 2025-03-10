/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Multiphase

  Self-registering factory for figacity models.
*/

#ifndef AMANZI_MULTIPHASE_FUGACITY_FACTORY_HH_
#define AMANZI_MULTIPHASE_FUGACITY_FACTORY_HH_

#include <memory>

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "Fugacity_HenryLaw.hh"
#include "Fugacity_IdealGas.hh"
#include "Fugacity_SaturatedVapor.hh"

namespace Amanzi {
namespace Multiphase {

class FugacityFactory : public Utils::Factory<Fugacity> {
 public:
  FugacityFactory(){};
  ~FugacityFactory(){};

  std::shared_ptr<Fugacity> Create(const Teuchos::ParameterList& plist) {
    std::string eos = plist.get<std::string>("eos type");
    if (eos == "ideal gas") {
      return std::make_shared<Fugacity_IdealGas>(plist);
    } else if (eos == "saturated vapor") {
      return std::make_shared<Fugacity_SaturatedVapor>(plist);
    } else if (eos == "Henry law") {
      return std::make_shared<Fugacity_SaturatedVapor>(plist);
    }
    return nullptr;
  }
};

} // namespace Multiphase
} // namespace Amanzi

#endif
