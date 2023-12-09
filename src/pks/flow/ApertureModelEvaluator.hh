/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  The fracture aperture evaluator.
*/

#ifndef AMANZI_FLOW_APERTURE_EVALUATOR_HH_
#define AMANZI_FLOW_APERTURE_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Key.hh"
#include "Tag.hh"

#include "ApertureModel.hh"
#include "ApertureModelPartition.hh"

namespace Amanzi {
namespace Flow {

class ApertureModelEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit ApertureModelEvaluator(Teuchos::ParameterList& plist,
                                  Teuchos::RCP<ApertureModelPartition> ppm);
  ApertureModelEvaluator(const ApertureModelEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<ApertureModelPartition> apm_;
  Key pressure_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, ApertureModelEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
