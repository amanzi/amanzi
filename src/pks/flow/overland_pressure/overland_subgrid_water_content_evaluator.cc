/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*

  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! An evaluator for determining the water content of a given pressure in a cell that includes subgrid microtopography.
/*!

This is Phi * density * area, where Phi is the volumetric (effective) ponded depth as in Jan et al WRR 2018.

* `"maximum ponded depth key`" ``[string]`` **DOMAIN-maximum_ponded_depth**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"pressure key`" ``[string]`` **DOMAIN-pressure**
         The name of the pressure on the surface.
* `"cell volume key`" ``[string]`` **DOMAIN-cell_volume**
         The name of the cell's volume.

*/

#include "overland_subgrid_water_content_evaluator.hh"

namespace Amanzi {
namespace Flow {

// Constructor from ParameterList
OverlandSubgridWaterContentEvaluator::OverlandSubgridWaterContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist){ 

  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  M_ = plist_.get<double>("molar mass", 0.0180153);
  n_liq_ = plist_.get<double>("molar density", 1000./0.0180153);
  //bar_ = plist_.get<bool>("allow negative water content", false);

  //
  // NOTE: A bar in this evaluator doesn't make sense, because bar is intended
  // to keep the derivative positive even for negative pressures, but this
  // curve asympotes to 0 with a 0 derivative, so negative values should be
  // zero.  That said I'm afraid to delete the code in case I'm missing
  // something.  FIXME!  --etc
  bar_ = false;
  
  Key domain = Keys::getDomain(my_key_);

  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief"); 
  dependencies_.insert(delta_max_key_);
  
  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(delta_ex_key_);

  pres_key_ = Keys::readKey(plist_, domain, "pressure", "pressure");
  dependencies_.insert(pres_key_);

  cv_key_ = Keys::readKey(plist_, domain, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
}


void
OverlandSubgridWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto& res = *result->ViewComponent("cell",false);

  const auto& pres = *S->GetFieldData(pres_key_)->ViewComponent("cell",false);
  const auto& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell",false);

  // these dependencies are assumed here -- must be updated in state-dev --etc
  const auto& p_atm = *S->GetScalarData("atmospheric_pressure");
  const auto& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];
  
  const auto& max_pd_v = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const auto& ex_vol_v = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);

  // if (bar_) {
  //   for (int c=0; c!=res.MyLength(); ++c) {
  //     // convert these parameters to water content
  //     double delta_max = max_pd_v[0][c] * n_liq_;
  //     double delta_ex = ex_vol_v[0][c] * n_liq_;

  //     double pd = (pres[0][c] - p_atm) / (gz * M_);
  //     if (pd <=delta_max){
  //       double pd_on_delmax =  pd / delta_max;
  //       res[0][c] = std::pow(pd_on_delmax,2) * (2*delta_max - 3*delta_ex)
  //                   + std::pow(pd_on_delmax,3) * (2*delta_ex - delta_max);
  //       res[0][c] *= cv[0][c];
  //     } else if (pd > delta_max) {
  //       res[0][c] = cv[0][c] *(pd - delta_ex);
  //     }
  //   }
  // } else {
    for (int c=0; c!=res.MyLength(); ++c) {
      double delta_max = max_pd_v[0][c] * n_liq_;
      double delta_ex = ex_vol_v[0][c] * n_liq_;

      double pd = (pres[0][c] - p_atm)/ (gz * M_);
      if (pd <= 0.) {
        res[0][c] = 0.;
      } else if (pd < delta_max) {
        double pd_on_delmax = pd / delta_max;
        res[0][c] = std::pow(pd_on_delmax,2) * (2*delta_max - 3*delta_ex)
                    + std::pow(pd_on_delmax,3) * (2*delta_ex - delta_max);
	res[0][c] *= cv[0][c];
      } else {
	res[0][c] = cv[0][c] * (pd - delta_ex);
      }
    }
  // }
}


void
OverlandSubgridWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(wrt_key == pres_key_);

  auto& res = *result->ViewComponent("cell",false);
  const auto& pres = *S->GetFieldData(pres_key_)->ViewComponent("cell",false);
  const auto& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell",false);

  const auto& p_atm = *S->GetScalarData("atmospheric_pressure");
  const auto& gravity = *S->GetConstantVectorData("gravity");
  double gz = -gravity[2];
 
  const auto& max_pd_v = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const auto& ex_vol_v = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);

  if (wrt_key == pres_key_) {
    // if (bar_) {
    //   for (int c=0; c!=res.MyLength(); ++c) {
    //     // This doesn't keep the value from 0, and I don't think it should have to in normal usage. --etc
    //     double pd = pres[0][c] < p_atm ? 1.0e-5 : (pres[0][c] - p_atm)/ (gz * M_);
    //     //double pd = (pres[0][c] - p_atm)/ (gz * M_);

    //     double delta_max = max_pd_v[0][c] * n_liq_;
    //     double delta_ex = ex_vol_v[0][c] * n_liq_;
    //     if (pd < delta_max){
    //       res[0][c] = 2 * pd * (2*delta_max - 3*delta_ex) / std::pow(delta_max,2)
    //                   + 3 * std::pow(pd, 2) * (2*delta_ex - delta_max) / std::pow(delta_max,3);
    //       res[0][c] *= cv[0][c] / (gz * M_);
    //     } else {
    //       res[0][c] = cv[0][c] / (gz * M_);
    //     }
    //   }

    // } else {
      for (int c=0; c!=res.MyLength(); ++c) {
        double pd = (pres[0][c] - p_atm) / (gz * M_);
        double delta_max =  max_pd_v[0][c] * n_liq_;
        double delta_ex = ex_vol_v[0][c] * n_liq_;
        if (pd <= 0.) {
          res[0][c] = 0.;
        } else if (pd < delta_max) {
          res[0][c] = 2 * pd * (2*delta_max - 3*delta_ex) / std::pow(delta_max,2)
                      + 3 * std::pow(pd, 2) * (2*delta_ex - delta_max) / std::pow(delta_max,3);
          res[0][c] *= cv[0][c] / (gz * M_);
        } else {
          res[0][c] = cv[0][c] / (gz * M_);
        }
      }
    // }
  } else {
    res.PutScalar(0.);
  }
}


} //namespace
} //namespace

