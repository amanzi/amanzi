/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Calculates liquid mole fraction from partial gas pressure.
*/

#ifndef AMANZI_MULTIPHASE_MOLAR_FRACTION_LIQUID_HH_
#define AMANZI_MULTIPHASE_MOLAR_FRACTION_LIQUID_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi

// Multiphase
#include "MultiphaseBaseEvaluator.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class TccLiquid : public MultiphaseBaseEvaluator {
 public:
  TccLiquid(Teuchos::ParameterList& plist);

  // interface methods from FieldEvaluator
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  std::string tcc_gas_key_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif

