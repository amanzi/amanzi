/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  MultiPhase PK

  Calculates mass density of gas phase.
*/

#ifndef AMANZI_MULTIPHASE_MASS_DENSITY_GAS_HH_
#define AMANZI_MULTIPHASE_MASS_DENSITY_GAS_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"

// Multiphase
#include "MultiphaseEvaluator.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class MassDensityGas : public MultiphaseEvaluator {
 public:
  MassDensityGas(Teuchos::ParameterList& plist);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  Key x_vapor_key_, x_gas_key_, mol_density_gas_key_;
  double mol_mass_H2O_;
  std::vector<double> mol_mass_;

  static Utils::RegisteredFactory<Evaluator, MassDensityGas> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
