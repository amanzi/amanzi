/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  Evaluator for computing shear modulus. Currently, the only available
  option is the Hardin-Drnevich model.
*/

#ifndef AMANZI_MECHANICS_SHEAR_MODULUS_EVALUATOR_
#define AMANZI_MECHANICS_SHEAR_MODULUS_EVALUATOR_

#include "PDE_Elasticity.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Mechanics {

class ShearStrainEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  ShearStrainEvaluator(Teuchos::ParameterList& plist);
  ShearStrainEvaluator(const ShearStrainEvaluator& other);

  // special initialization function
  void set_op(Teuchos::RCP<Operators::PDE_Elasticity>& op) { op_ = op; }

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

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
};

} // namespace Mechanics
} // namespace Amanzi

#endif
