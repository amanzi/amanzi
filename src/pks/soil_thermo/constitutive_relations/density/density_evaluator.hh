/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for water density.
----------------------------------------------------------------------------- */


#ifndef AMANZI_SOIL_DENSITY_EVALUATOR_HH_
#define AMANZI_SOIL_DENSITY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

class DensityEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  DensityEvaluator(Teuchos::ParameterList& plist);
  DensityEvaluator(const DensityEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:

  Key temperature_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,DensityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
