/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for energy, e = cp*rho*T.
----------------------------------------------------------------------------- */


#ifndef AMANZI_SOIL_ENERGY_EVALUATOR_HH_
#define AMANZI_SOIL_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

class SoilEnergyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SoilEnergyEvaluator(Teuchos::ParameterList& plist);
  SoilEnergyEvaluator(const SoilEnergyEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key temperature_key_;
  Key pres_key_;
  Key density_key_;
  Key heat_capacity_key_;
  Key pressure_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SoilEnergyEvaluator> factory_;

};

} // namespace
} // namespace

#endif
