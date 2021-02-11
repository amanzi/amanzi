/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "boost/algorithm/string/predicate.hpp"

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
  min_area_ = plist_.get<double>("minimum fractional area [-]", 1.e-5);
  if (min_area_ <= 0.) {
    Errors::Message message("AreaFractionsEvaluator: Minimum fractional area should be > 0.");
    Exceptions::amanzi_throw(message);
  }

  // get domain names
  domain_ = Keys::getDomain(my_key_);
  if (domain_ == "surface") {
    domain_snow_ = plist_.get<std::string>("snow domain name", "snow");
  } else if (boost::starts_with(domain_, "surface_")) {
    domain_snow_ = plist_.get<std::string>("snow domain name",
            std::string("snow_") + domain_.substr(8,domain_.size()));
  } else {
    domain_snow_ = plist_.get<std::string>("snow domain name");
  }

  // get dependencies
  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(snow_depth_key_);
}


void
AreaFractionsEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto mesh = S->GetMesh(domain_);

  auto& res = *result->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);

  for (const auto& lc : land_cover_) {
    Entity_ID_List lc_ids;
    mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_type::OWNED, lc_ids);

    for (auto c : lc_ids) {
      // calculate area of land
      if (sd[0][c] >= lc.second.snow_transition_depth) {
        res[0][c] = 0.;
        res[1][c] = 1.;
      } else if (sd[0][c] <= 0.) {
        res[0][c] = 1.;
        res[1][c] = 0.;
      } else {
        res[1][c] = sd[0][c] / lc.second.snow_transition_depth;
        res[0][c] = 1. - res[1][c];
      }

      // if any area is less than eps, give to other
      if (res[0][c] > 0 && res[0][c] < min_area_) {
        res[0][c] = 0.;
        res[1][c] = 1.;
      } else if (res[1][c] > 0 && res[1][c] < min_area_) {
        res[0][c] = 1.;
        res[1][c] = 0.;
      }

      AMANZI_ASSERT(std::abs(res[0][c] + res[1][c] - 1.0) < 1.e-10);
      AMANZI_ASSERT(-1.e-10 <= res[0][c] && res[0][c] <= 1.+1.e-10);
      AMANZI_ASSERT(-1.e-10 <= res[1][c] && res[1][c] <= 1.+1.e-10);

      res[0][c] = std::min(std::max(0.,res[0][c]), 1.);
      res[1][c] = std::min(std::max(0.,res[1][c]), 1.);
    }
  }
}

void
AreaFractionsEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // new state!
  land_cover_ = getLandCover(S->ICList()->sublist("land cover types"));

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
    if (Keys::getDomain(dep_key) == domain_snow_) {
      fac->SetMesh(S->GetMesh(domain_snow_))
          ->SetGhosted()
          ->SetComponent("cell", AmanziMesh::CELL, 1);
    } else {
      fac->SetMesh(S->GetMesh(domain_))
          ->SetGhosted()
          ->SetComponent("cell", AmanziMesh::CELL, 1);
    }

    // Recurse into the tree to propagate info to leaves.
    S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
  }
}


} //namespace
} //namespace

