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

* `"snow depth key`" ``[string]`` **DOMAIN-depth**
         The name of the snow depth.

* `"snow-ground transitional depth [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the snow gradient.
         
Ordering of the area fractions calculated are: [land, snow].
         
NOTE: this evaluator simplifies the situation by assuming constant
density.  This make it so that ice and water see the same geometry per
unit pressure, which isn't quite true thanks to density differences.
However, we hypothesize that these differences, on the surface (unlike in
the subsurface) really don't matter much. --etc
         
*/

#include "area_fractions_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// Constructor from ParameterList
AreaFractionsEvaluator::AreaFractionsEvaluator(Teuchos::ParameterList& plist) :
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
  if (domain_ == "surface") {
    domain_snow_ = "snow";
  } else if (boost::starts_with(domain_, "surface_")) {
    domain_snow_ = std::string("snow_") + domain_.substr(8,domain_.size());
  }

  // FIXME: "maximum_ponded_depth" is a terrible name, this is a geometric thing, not a dynamic thing. --etc
  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(snow_depth_key_);
}


void
AreaFractionsEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto& res = *result->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);

  for (int c=0; c!=res.MyLength(); ++c) {
    // calculate area of land
    if (sd[0][c] >= snow_subgrid_transition_) {
      res[0][c] = 0.;
      res[1][c] = 1.;
    } else if (sd[0][c] <= min_area_) {
      res[0][c] = 1.;
      res[1][c] = 0.;
    } else {
      res[1][c] = (sd[0][c] - min_area_) / (snow_subgrid_transition_ - min_area_);
      res[0][c] = 1. - res[1][c];
    }
  }
}

void
AreaFractionsEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // see if we can find a master fac
  auto my_fac = S->RequireField(my_key_, my_key_);
  my_fac->SetMesh(S->GetMesh(domain_))
      ->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 2);

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

