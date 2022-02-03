/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Determining the molar fraction of a gas component within a gas mixture.
*/

#ifndef AMANZI_EOS_MOLAR_FRACTION_GAS_EVALUATOR_HH_
#define AMANZI_EOS_MOLAR_FRACTION_GAS_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

#include "EOS_SaturatedVaporPressure.hh"

namespace Amanzi {
namespace AmanziEOS {

class MolarFractionGasEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit MolarFractionGasEvaluator(Teuchos::ParameterList& plist);
  MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  Teuchos::RCP<EOS_SaturatedVaporPressure> get_svp_model() { return svp_model_; }

 protected:
  Key temp_key_;
  Tag tag_;

  Teuchos::RCP<EOS_SaturatedVaporPressure> svp_model_;

 private:
  static Utils::RegisteredFactory<Evaluator, MolarFractionGasEvaluator> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
