/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a soil heat capacity model

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SOIL_HC_EVALUATOR_HH_
#define AMANZI_SOIL_HC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

class SoilHeatCapacityEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  SoilHeatCapacityEvaluator(Teuchos::ParameterList& plist);
  SoilHeatCapacityEvaluator(const SoilHeatCapacityEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // dependencies

  double cg, cw, ci;

  Key temperature_key_;
  Key water_content_key_;
  Key ice_content_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SoilHeatCapacityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
