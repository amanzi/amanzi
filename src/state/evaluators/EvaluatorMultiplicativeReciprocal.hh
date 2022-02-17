/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL. 
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

#include <string>

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class EvaluatorMultiplicativeReciprocal : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  EvaluatorMultiplicativeReciprocal(Teuchos::ParameterList& plist);
  EvaluatorMultiplicativeReciprocal(const EvaluatorMultiplicativeReciprocal& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 private:
  bool enforce_positivity_;
  double coef_;
  std::vector<std::string> list0_, list1_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorMultiplicativeReciprocal> fac_;
};

}  // namespace Amanzi

#endif
