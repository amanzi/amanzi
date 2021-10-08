/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Secondary variable field evaluator computes product of fields 
  or inverse of fields:

    eval = f1 * f2 * ... * fn) / (g1 * g2 * ... * gm)
*/

#ifndef AMANZI_MULTIPLICATIVE_RECIPROCAL_EVALUATOR_HH_
#define AMANZI_MULTIPLICATIVE_RECIPROCAL_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

class MultiplicativeReciprocalEvaluator : public SecondaryVariableFieldEvaluator {
 public:
  MultiplicativeReciprocalEvaluator(Teuchos::ParameterList& plist);
  MultiplicativeReciprocalEvaluator(const MultiplicativeReciprocalEvaluator& other);

  // inteface functions to FieldEvaluator
  virtual Teuchos::RCP<FieldEvaluator> Clone() const override;

  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
      Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

 private:
  std::vector<std::string> list0_, list1_;
};

}  // namespace Amanzi

#endif
