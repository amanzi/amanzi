/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "surface_coupler_via_source_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

Utils::RegisteredFactory<FieldEvaluator,SurfaceCouplerViaSourceEvaluator>
SurfaceCouplerViaSourceEvaluator::fac_("surface_coupler_via_source");

} // namespace
} // namespace
} // namespace
