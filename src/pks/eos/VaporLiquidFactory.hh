/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Factory for vapor-liquid distribution coefficient.
*/

#ifndef AMANZI_EOS_VAPOR_LIQUID_FACTORY_HH_
#define AMANZI_EOS_VAPOR_LIQUID_FACTORY_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"

#include "VaporLiquid_Constant.hh"
#include "VaporLiquid_Tabular.hh"

namespace Amanzi {
namespace AmanziEOS {

class VaporLiquidFactory {
 public:
  VaporLiquidFactory(const Teuchos::ParameterList& plist) : plist_(plist) {};

  std::shared_ptr<VaporLiquid> Create(const std::string& name) {
    if (EquilibriumCoef.find(name) != EquilibriumCoef.end()) {
      return std::make_shared<VaporLiquid_Tabular>(name);
    } else {
      const auto& names = plist_.get<Teuchos::Array<std::string>>("aqueous names").toVector();
      const auto& values = plist_.get<Teuchos::Array<double>>("Henry dimensionless constants").toVector();
     
      for (int i = 0; i < names.size(); ++i) {
        if (names[i] == name)
          return std::make_shared<VaporLiquid_Constant>(values[i]);
      }
    }

    AMANZI_ASSERT(false);
    return nullptr;
  }

 private:
  Teuchos::ParameterList plist_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
