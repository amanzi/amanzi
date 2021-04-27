/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of land and snow within a grid cell

/*!

Requires the following dependencies:

* `"snow depth key`" ``[string]`` **DOMAIN-snow_depth**

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

class AreaFractionsThreeChannelEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  AreaFractionsThreeChannelEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsThreeChannelEvaluator(const AreaFractionsThreeChannelEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new AreaFractionsThreeChannelEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Exceptions::amanzi_throw("NotImplemented: AreaFractionsThreeChannelEvaluator currently does not provide derivatives.");
  }

 protected:

  Key domain_, domain_snow_;
  Key snow_depth_key_;
  Key ponded_depth_key_;
  double snow_subgrid_transition_;
  double water_subgrid_transition_;
  double min_area_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AreaFractionsThreeChannelEvaluator> reg_;

};

} //namespace
} //namespace


