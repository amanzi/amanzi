/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! An evaluator with no dependencies specified by a function of t,x,y,z.

/*!

.. todo:
    This needs a test and documentation! --etc

*/


#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFUNCTION_

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentFunction
  : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentFunction(Teuchos::ParameterList& plist);
  EvaluatorIndependentFunction(const EvaluatorIndependentFunction& other) =
    default;

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
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction> fac_;
};

} // namespace Amanzi

#endif
