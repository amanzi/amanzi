/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License: see COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies specified by a function.

TODO: This needs a test! --etc
------------------------------------------------------------------------- */

#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentFunction :
      public EvaluatorIndependent<CompositeVector,CompositeVectorSpace> {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  EvaluatorIndependentFunction(Teuchos::ParameterList& plist);
  EvaluatorIndependentFunction(const EvaluatorIndependentFunction& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentFunction&
  operator=(const EvaluatorIndependentFunction& other);

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> func_;

 private:
  static Utils::RegisteredFactory<Evaluator,EvaluatorIndependentFunction> fac_;
};

} // namespace


#endif
