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

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#ifndef AMANZI_FLOW_WRM_EVALUATOR_HH_
#define AMANZI_FLOW_WRM_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
#include "Tag.hh"

#include "WRM.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

class WRMEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit WRMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrm);
  WRMEvaluator(const WRMEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  void InitializeFromPlist_();

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<WRMPartition> wrm_;
  Key pressure_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, WRMEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
