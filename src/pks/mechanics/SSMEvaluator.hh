/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  The evaluator calls the small strain model with the correct arguments.
*/

#ifndef AMANZI_MECHANICS_SSM_EVALUATOR_HH_
#define AMANZI_MECHANICS_SSM_EVALUATOR_HH_

#include <utility>

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Tag.hh"

#include "SSM.hh"
#include "SSMPartition.hh"

namespace Amanzi {
namespace Mechanics {

class SSMEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit SSMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<SSMPartition>& ssm);
  SSMEvaluator(const SSMEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<SSMPartition> ssm_;
  Key shear_strain_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SSMEvaluator> reg_;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
