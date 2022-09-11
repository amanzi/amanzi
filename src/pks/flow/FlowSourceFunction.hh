/*
  Flow PK
 
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_SOURCE_FUNCTION_HH_
#define AMANZI_FLOW_SOURCE_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "PK_DomainFunction.hh"
#include "State.hh"

namespace Amanzi {
namespace Flow {

class FlowSourceFunction : public PK_DomainFunction {
 public:
  FlowSourceFunction() {};
  FlowSourceFunction(const Teuchos::ParameterList& plist) {};

  void ComputeSubmodel(const Key& key, const State& S) {
    if (name() != "volume" && S.HasRecord(key, Tags::DEFAULT)) {
      auto aperture = *S.Get<CompositeVector>(key, Tags::DEFAULT).ViewComponent("cell", true);
      for (auto it = begin(); it != end(); ++it) {
        it->second[0] *= aperture[0][it->first];
      }
    }
  }
};

}  // namespace Flow
}  // namespace Amanzi

#endif
