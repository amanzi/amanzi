/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "surface_relperm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SurfaceRelPermEvaluator>
SurfaceRelPermEvaluator::fac_("surface rel perm");

} //namespace
} //namespace
} //namespace

