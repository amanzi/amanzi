/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A subgrid model for determining the area fraction of land and snow within a grid cell
/*!

A simple linear scaling from bare to snow-covered area fractions [-].  Ordering
of the area fractions calculated are: [land, snow].

.. _area-fractions-evaluator-spec:
.. admonition:: area-fractions-evaluator-spec

   * `"snow transitional height [m]`" ``[double]`` **0.02**
     Minimum thickness for specifying the snow gradient.
   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
     Mimimum area fraction allowed, less than this is rebalanced as zero.
   * `"snow domain name`" ``[string]`` **DOMAIN_SNOW** A default is guessed at
     by replacing `"surface`" with `"snow`" in the this's domain.

   KEYS:
   * `"snow depth`" **DOMAIN_SNOW-depth**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

class AreaFractionsEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  AreaFractionsEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsEvaluator(const AreaFractionsEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new AreaFractionsEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Exceptions::amanzi_throw("NotImplemented: AreaFractionsEvaluator currently does not provide derivatives.");
  }

 protected:

  Key domain_, domain_snow_;
  Key snow_depth_key_;
  double min_area_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverVector land_cover_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AreaFractionsEvaluator> reg_;

};

} //namespace
} //namespace
