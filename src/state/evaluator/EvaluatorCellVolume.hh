/*
  State

  License: BSD
  Author: Ethan Coon

*/

#ifndef AMANZI_INDEPENDENT_CELL_VOLUME_EVALUATOR_HH_
#define AMANZI_INDEPENDENT_CELL_VOLUME_EVALUATOR_HH_

#include "EvaluatorIndependent.hh"

namespace Amanzi {

//
// Dummy class, does everything but know the type, which is required to
// EnsureCompatibility.  This is never used, instead the below templated one
// is.
//
class EvaluatorCellVolume : public EvaluatorIndependent<CompositeVector,CompositeVectorSpace> {

public:
  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  using EvaluatorIndependent<CompositeVector,CompositeVectorSpace>::EvaluatorIndependent;
  EvaluatorCellVolume(const EvaluatorCellVolume &other) = default;

  Teuchos::RCP<Evaluator> Clone() const {
    return Teuchos::rcp(new EvaluatorCellVolume(*this));
  }

protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State &S);

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorCellVolume> fac_;
};

} // namespace Amanzi

#endif
