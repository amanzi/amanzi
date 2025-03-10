/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  Determining the molar fraction of a gas component within a gas mixture.
*/

#ifndef AMANZI_ENERGY_MOLE_FRACTION_GAS_EVALUATOR_HH_
#define AMANZI_ENERGY_MOLE_FRACTION_GAS_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

#include "EOS_SaturatedVaporPressure.hh"

namespace Amanzi {
namespace Energy {

class MoleFractionGasEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit MoleFractionGasEvaluator(Teuchos::ParameterList& plist);
  MoleFractionGasEvaluator(const MoleFractionGasEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Key temp_key_;
  Tag tag_;

  Teuchos::RCP<AmanziEOS::EOS_SaturatedVaporPressure> svp_model_;

 private:
  static Utils::RegisteredFactory<Evaluator, MoleFractionGasEvaluator> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
