/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of land, water, and snow within a grid cell with subgrid microtopography.

/*!

  Uses the subgrid equation from Jan et al WRR 2018 for volumetric or
  effective ponded depth to determine the area of water, then heuristically
  places snow on top of that surface.

Requires the following dependencies:

* `"microtopographic relief key`" ``[string]`` **DOMAIN-maximum_ponded_depth**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-pressure**
         The name of the surface water ponded depth.
* `"snow depth key`" ``[string]`` **DOMAIN-snow_depth**
         The name of the snow depth.

* `"snow-ground transitional depth [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the snow gradient.
         
Ordering of the area fractions calculated are: [land, water, snow].
         
NOTE: this evaluator simplifies the situation by assuming constant
density.  This make it so that ice and water see the same geometry per
unit pressure, which isn't quite true thanks to density differences.
However, we hypothesize that these differences, on the surface (unlike in
the subsurface) really don't matter much. --etc
         
*/

#include "area_fractions_subgrid_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// Constructor from ParameterList
AreaFractionsSubgridEvaluator::AreaFractionsSubgridEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{ 
  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  snow_subgrid_transition_ = plist_.get<double>("snow transition height [m]", 0.02);
  min_area_ = plist_.get<double>("minimum fractional area [-]", 0.);
  
  domain_ = Keys::getDomain(my_key_);

  // FIXME: "maximum_ponded_depth" is a terrible name, this is a geometric thing, not a dynamic thing. --etc
  delta_max_key_ = Keys::readKey(plist_, domain_, "microtopographic relief", "maximum_ponded_depth"); 
  dependencies_.insert(delta_max_key_);
  
  delta_ex_key_ = Keys::readKey(plist_, domain_, "excluded volume", "excluded_volume");
  dependencies_.insert(delta_ex_key_);

  ponded_depth_key_ = Keys::readKey(plist_, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(ponded_depth_key_);

  snow_depth_key_ = Keys::readKey(plist_, domain_, "snow depth", "snow_depth");
  dependencies_.insert(snow_depth_key_);

  vol_snow_depth_key_ = Keys::readKey(plist_, domain_, "volumetric snow depth", "volumetric_snow_depth");
  dependencies_.insert(vol_snow_depth_key_);
}


void
AreaFractionsSubgridEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto& res = *result->ViewComponent("cell",false);

  const auto& pd = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);
  const auto& vsd = *S->GetFieldData(vol_snow_depth_key_)->ViewComponent("cell",false);
  const auto& del_max = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const auto& del_ex = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);

  for (int c=0; c!=res.MyLength(); ++c) {
    // calculate area of land
    double wet_area = f_prime_(pd[0][c] + sd[0][c], del_max[0][c], del_ex[0][c]);
    res[0][c] = 1 - wet_area;

    // now partition the wet area into snow and water
    if (vsd[0][c] >= wet_area * snow_subgrid_transition_) {
      res[2][c] = wet_area;
      res[1][c] = 0.;
    } else {
      res[2][c] = vsd[0][c] / snow_subgrid_transition_;
      res[1][c] = wet_area - res[2][c];
    }

    AMANZI_ASSERT(std::abs(res[0][c] + res[1][c] + res[2][c] - 1.0) < 1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[0][c] && res[0][c] <= 1.+1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[1][c] && res[1][c] <= 1.+1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[2][c] && res[2][c] <= 1.+1.e-10);
    res[0][c] = res[0][c] < min_area_ ? 0. : std::min(std::max(0.,res[0][c]), 1.);
    res[1][c] = res[1][c] < min_area_ ? 0. : std::min(std::max(0.,res[1][c]), 1.);
    res[2][c] = res[2][c] < min_area_ ? 0. : std::min(std::max(0.,res[2][c]), 1.);
  }
}

void
AreaFractionsSubgridEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // see if we can find a master fac
  auto my_fac = S->RequireField(my_key_, my_key_);
  my_fac->SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Check plist for vis or checkpointing control.
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  
  for (auto dep_key : dependencies_) {
    auto fac = S->RequireField(dep_key);
    fac->SetMesh(S->GetMesh(domain_))
        ->SetGhosted()
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    // Recurse into the tree to propagate info to leaves.
    S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
  }
}

} //namespace
} //namespace

