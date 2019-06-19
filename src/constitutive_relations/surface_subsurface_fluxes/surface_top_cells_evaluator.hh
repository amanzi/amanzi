/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_SURFACE_TOP_CELLS_EVALUATOR_
#define AMANZI_RELATIONS_SURFACE_TOP_CELLS_EVALUATOR_

#include "Factory.hh"

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class SurfaceTopCellsEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SurfaceTopCellsEvaluator(Teuchos::ParameterList& plist);

  SurfaceTopCellsEvaluator(const SurfaceTopCellsEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    AMANZI_ASSERT(0);
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  Key dependency_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SurfaceTopCellsEvaluator> reg_;

};

} //namespace
} //namespace

#endif
