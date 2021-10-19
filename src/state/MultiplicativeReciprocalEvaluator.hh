/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  Secondary variable field evaluator computes product of fields 
  or inverse of fields:

    eval = f1 * f2 * ... * fn) / (g1 * g2 * ... * gm)
*/

//! A generic evaluator for multiplying and dividing a collection of fields.

/*!

.. _multiplicative-reciprocal-evaluator-spec:
.. admonition:: multiplicative-reciprocal-evaluator-spec
   * `"coefficient`" ``[double]`` **1** A constant prefix to the product.
   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   * `"multiplicative dependencies`" ``[Array(string)]`` **optional**, only base field names
   * `"reciprocal dependencies`" ``[Array(string)]`` **optional**, only base field names

*/
#ifndef AMANZI_MULTIPLICATIVE_RECIPROCAL_EVALUATOR_HH_
#define AMANZI_MULTIPLICATIVE_RECIPROCAL_EVALUATOR_HH_

#include "Factory.hh"
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
  bool enforce_positivity_;
  double coef_;
  std::vector<std::string> list0_, list1_;

  static Utils::RegisteredFactory<FieldEvaluator, MultiplicativeReciprocalEvaluator> factory_;
};

}  // namespace Amanzi

#endif
