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
#include "Factory.hh"

// Multiphase
#include "MultiphaseBaseEvaluator.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class MoleFractionLiquid : public MultiphaseBaseEvaluator {
 public:
  MoleFractionLiquid(Teuchos::ParameterList& plist);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  std::string x_gas_key_, pressure_gas_key_;

  static Utils::RegisteredFactory<Evaluator, MoleFractionLiquid> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
