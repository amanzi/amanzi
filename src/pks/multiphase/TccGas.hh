/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Calculates gas component concentration from liquid concentration.
*/

#ifndef AMANZI_MULTIPHASE_TCC_GAS_HH_
#define AMANZI_MULTIPHASE_TCC_GAS_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"

// Multiphase
#include "MultiphaseBaseEvaluator.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class TccGas : public MultiphaseBaseEvaluator {
 public:
  TccGas(Teuchos::ParameterList& plist);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  std::string tcc_liquid_key_;

  static Utils::RegisteredFactory<Evaluator, TccGas> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
