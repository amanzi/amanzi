/*
  Shallow Water PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstanrtin Lipnikov

  Self-registering factory for numerical flux implementations.
*/

#ifndef SHALLOW_WATER_NUMERICAL_FLUX_FACTORY_HH_
#define SHALLOW_WATER_NUMERICAL_FLUX_FACTORY_HH_

#include <memory>

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "NumericalFlux.hh"

namespace Amanzi {
namespace ShallowWater {

class NumericalFluxFactory : public Utils::Factory<NumericalFlux> {
 public:
  std::shared_ptr<NumericalFlux> Create(Teuchos::ParameterList& plist)
  {
    std::string model = plist.get<std::string>("numerical flux");
    return std::shared_ptr<NumericalFlux>(CreateInstance(model, plist));
  }
};

} // namespace ShallowWater
} // namespace Amanzi

#endif
