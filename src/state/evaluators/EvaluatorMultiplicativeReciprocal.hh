/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
/*!

A generic evaluator for multiplying and dividing a collection of fields.

Secondary variable field evaluator computes element-wise multiplication and
division of fields:

.. math::
   e = \alpha \frac{f_1 * f_2 * ... * f_n}{g_1 * g_2 * ... * g_m}

Note that, for the moment, all tags for :math:`f_i` and :math:`g_i` must be the
same as that of :math:`e`.  This could straightforwardly be relaxed upon
request.

`"evaluator type`" = `"multiplicative reciprocal`"

.. _evaluator-multiplicative-reciprocal-spec:
.. admonition:: evaluator-multiplicative-reciprocal-spec

   * `"coefficient`" ``[double]`` **1** A constant prefix to the product.
   * `"enforce positivity`" ``[bool]`` **false** If true, max the result with 0.

   ONE OF

   * `"multiplicative dependency keys`" ``[Array(string)]`` **optional**, the :math:`f_i`

   OR

   * `"multiplicative dependency key suffixes`" ``[Array(string)]`` **optional**, suffixes for the :math:`f_i`.

   END

   ONE OF

   * `"reciprocal dependency keys`" ``[Array(string)]`` **optional**, the :math:`g_i`

   OR

   * `"reciprocal dependency key suffixes`" ``[Array(string)]`` **optional**, suffixes for the :math:`g_i`.

   END

*/

#pragma once

#include <string>

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class EvaluatorMultiplicativeReciprocal
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  EvaluatorMultiplicativeReciprocal(Teuchos::ParameterList& plist);
  EvaluatorMultiplicativeReciprocal(const EvaluatorMultiplicativeReciprocal& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_Units_(State& S) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 private:
  bool enforce_positivity_;
  double coef_;
  Teuchos::Array<std::string> list0_, list1_;
  int n_dofs_;

  static Utils::RegisteredFactory<Evaluator, EvaluatorMultiplicativeReciprocal> fac_;
};

} // namespace Amanzi
