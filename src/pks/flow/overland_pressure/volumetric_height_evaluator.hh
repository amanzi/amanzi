/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_VOLUMETRIC_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_VOLUMETRIC_HEIGHT_EVALUATOR_

//#include "height_evaluator.hh"
#include "secondary_variable_field_evaluator.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {

class VolumetricHeightModel;

class VolumetricHeightEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  VolumetricHeightEvaluator(Teuchos::ParameterList& plist);
  VolumetricHeightEvaluator(const VolumetricHeightEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  Teuchos::RCP<VolumetricHeightModel> get_VolumetricModel() { return vol_model_; }

 protected:
//  void InitializeFromPlist_();

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key pd_key_;
  Key delta_max_key_, delta_ex_key_;
  Teuchos::RCP<VolumetricHeightModel> vol_model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,VolumetricHeightEvaluator> factory_;

};

} //namespace
} //namespace

#endif
