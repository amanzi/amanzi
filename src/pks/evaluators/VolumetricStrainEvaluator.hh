/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

  Evaluator for computing hydrostatic stress.
*/

#ifndef AMANZI_VOLUMETRIC_STRAIN_EVALUATOR_
#define AMANZI_VOLUMETRIC_STRAIN_EVALUATOR_

#include "Factory.hh"
#include "PDE_Elasticity.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Evaluators {

class VolumetricStrainEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  VolumetricStrainEvaluator(Teuchos::ParameterList& plist);
  VolumetricStrainEvaluator(const VolumetricStrainEvaluator& other);

  // special initialization function
  void set_op(Teuchos::RCP<Operators::PDE_Elasticity>& op) { op_ = op; }

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void EnsureEvaluators(State& S) override {};

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override {};

  // since cell-based stress requires nodes, the default behavior is not applicable
  virtual void EnsureCompatibility_ToDeps_(State& S) override {};

 protected:
  Key displacement_key_;
  Teuchos::RCP<Operators::PDE_Elasticity> op_;

 private:
  static Utils::RegisteredFactory<Evaluator, VolumetricStrainEvaluator> reg_;
};

} // namespace Evaluators
} // namespace Amanzi

#endif
