/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

Interpolates between two time tags to compute a value at an intermediate time tag.

`"evaluator type`" == `"temporal interpolation`"

.. _temporal-interpolation-evaluator-spec:
.. admonition:: temporal-interpolation-evaluator-spec

   * `"current tag`" ``[string]`` Tag for the old or current time.
   * `"next tag`" ``[string]`` Tag for the new or next time.


This example is a case where transport needs to interpolate the water content
between the default, global CURRENT and NEXT tags to compute a (subcycled)
water content at an intermediate time (here called transport_current).

.. code-block:: xml

  <ParameterList name="water_content@transport_current">
    <Parameter name="current tag" type="string" value="CURRENT"/>
    <Parameter name="next tag" type="string" value=""/>
  </Parameter>

 */


#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class EvaluatorTemporalInterpolation : public EvaluatorSecondary {
 public:
  explicit EvaluatorTemporalInterpolation(Teuchos::ParameterList& plist);
  EvaluatorTemporalInterpolation(const EvaluatorTemporalInterpolation& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void EnsureCompatibility(State& S) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Update_(State& S) override;

  // could implement this, but harder because chain rule is not in
  // EvaluatorSecondary.  Do we need it?
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {}

  virtual bool IsDifferentiableWRT(const State& S,
                                   const Key& wrt_key,
                                   const Tag& wrt_tag) const override
  {
    return false;
  }

 protected:
  KeyTag current_, next_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorTemporalInterpolation> fac_;
};

} // namespace Amanzi
