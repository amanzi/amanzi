/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! An evaluator for determining the water content of a given pressure in a cell that includes subgrid microtopography.
/*!

This is Phi * density * area, where Phi is the volumetric (effective) ponded depth as in Jan et al WRR 2018.

* `"microtopographic relief key`" ``[string]`` **DOMAIN-microtopographic_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"pressure key`" ``[string]`` **DOMAIN-pressure**
         The name of the pressure on the surface.
* `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**
         The name of the cell's volume.

*/

#ifndef AMANZI_FLOW_OVERLAND_SUBGRID_WATER_CONTENT_EVALUATOR_HH_
#define AMANZI_FLOW_OVERLAND_SUBGRID_WATER_CONTENT_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class OverlandSubgridWaterContentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  OverlandSubgridWaterContentEvaluator(Teuchos::ParameterList& plist);
  OverlandSubgridWaterContentEvaluator(const OverlandSubgridWaterContentEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new OverlandSubgridWaterContentEvaluator(*this));
  }    

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key pres_key_, cv_key_;
  double n_liq_, M_;
  Key delta_max_key_, delta_ex_key_;
  bool bar_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,OverlandSubgridWaterContentEvaluator> reg_;

};

} //namespace
} //namespace

#endif
