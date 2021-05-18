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
namespace Relations {

// Constructor from ParameterList
AreaFractionsTwoComponentEvaluator::AreaFractionsTwoComponentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  min_area_ = plist_.get<double>("minimum fractional area [-]", 1.e-5);
  if (min_area_ <= 0.) {
    Errors::Message message("AreaFractionsTwoComponentEvaluator: Minimum fractional area should be > 0.");
    Exceptions::amanzi_throw(message);
  }

  // get domain names
  Key domain = Keys::getDomain(my_key_);

  // get dependencies
  snow_depth_key_ = Keys::readKey(plist_, domain, "snow depth", "depth");
  dependencies_.insert(snow_depth_key_);
}


void
AreaFractionsTwoComponentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"));

  auto mesh = result->Mesh();
  auto& res = *result->ViewComponent("cell",false);
  const auto& sd = *S->GetFieldData(snow_depth_key_)->ViewComponent("cell",false);

  for (const auto& lc : land_cover_) {
    AmanziMesh::Entity_ID_List lc_ids;
    mesh->get_set_entities(lc.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (auto c : lc_ids) {
      // calculate area of land
      if (sd[0][c] >= lc.second.snow_transition_depth) {
        res[0][c] = 1.;
      } else if (sd[0][c] <= 0.) {
        res[0][c] = 0.;
      } else {
        res[0][c] = sd[0][c] / lc.second.snow_transition_depth;
      }

      // if any area is less than eps, give to other
      if (res[0][c] < min_area_) {
        res[0][c] = 0.;
      } else if (res[0][c] > (1-min_area_)) {
        res[0][c] = 1.;
      }
    }
  }
}

void
AreaFractionsTwoComponentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Exceptions::amanzi_throw("NotImplemented: AreaFractionsTwoComponentEvaluator currently does not provide derivatives.");
}

} //namespace
} //namespace
} //namespace

