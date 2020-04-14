/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! An evaluator with no dependencies specified by a constant value.

/*!

.. todo:
    This needs a test and documentation! --etc

*/


#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_CONSTANT_HH_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_CONSTANT_HH_

#include "EvaluatorIndependent.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorIndependentConstant
  : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependentConstant(Teuchos::ParameterList& plist);
  EvaluatorIndependentConstant(const EvaluatorIndependentConstant& other) =
    default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorIndependentConstant&
  operator=(const EvaluatorIndependentConstant& other);

  virtual std::string name() const override { return "independent variable constant value"; }
  
 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

 protected:
  double value_;
  bool computed_once_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorIndependentConstant> fac_;
};

} // namespace Amanzi

#endif
