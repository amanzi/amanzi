/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SOIL_TC_EVALUATOR_HH_
#define AMANZI_SOIL_TC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

class SoilThermalConductivityEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  SoilThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  SoilThermalConductivityEvaluator(const SoilThermalConductivityEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  // dependencies

  Key temperature_key_;
  Key water_content_key_;
  Key ice_content_key_;
  Key cell_is_ice_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SoilThermalConductivityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
