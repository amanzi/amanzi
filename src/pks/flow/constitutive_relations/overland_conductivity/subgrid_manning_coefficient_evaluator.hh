/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! SubgridManningCoefficientEvaluator: calculates a fractional conductance.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan
*/

/*!

*/

#ifndef AMANZI_FLOWRELATIONS_SUBGRID_MANNING_COEF_EVALUATOR_HH_
#define AMANZI_FLOWRELATIONS_SUBGRID_MANNING_COEF_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class SubgridManningCoefficientEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SubgridManningCoefficientEvaluator(Teuchos::ParameterList& plist);

  SubgridManningCoefficientEvaluator(const SubgridManningCoefficientEvaluator& other) = default;

  Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new SubgridManningCoefficientEvaluator(*this));
  }
  
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 private:
  Key mann_key_, beta_key_, frac_cond_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SubgridManningCoefficientEvaluator> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
