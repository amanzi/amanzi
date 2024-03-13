/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Multiphase PK

  Field evaluator for noninear complimentary problem, function G.
*/

#ifndef AMANZI_MULTIPHASE_NCP_MOLAR_DENSITIES_HH_
#define AMANZI_MULTIPHASE_NCP_MOLAR_DENSITIES_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiphaseEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

class NCP_MolarDensities : public MultiphaseEvaluator {
 public:
  NCP_MolarDensities(Teuchos::ParameterList& plist);
  NCP_MolarDensities(const NCP_MolarDensities& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key x_vapor_key_, tcc_gas_key_, mol_density_gas_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, NCP_MolarDensities> fac_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
