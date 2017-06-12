/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "independent_variable_field_evaluator.hh"

#include "upwinding.hh"
#include "upwind_potential_difference.hh"
#include "pres_elev_evaluator.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "unfrozen_effective_depth_evaluator.hh"
#include "unfrozen_fraction_model.hh"
#include "unfrozen_fraction_evaluator.hh"

#include "overland_pressure_water_content_evaluator.hh"
#include "icy_overland.hh"

namespace Amanzi {
namespace Flow {

void IcyOverlandFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // ensure that the overland conductivity uses the unfrozen ponded depth
  // -- set the height key to be eta * h, not just h, for the frozen case.
  ASSERT(plist_->isSublist("overland conductivity evaluator"));

  if (!plist_->sublist("overland conductivity evaluator").isParameter("height key"))
    plist_->sublist("overland conductivity evaluator").set("height key", getKey(domain_,"unfrozen_effective_depth"));
  ASSERT(plist_->sublist("overland conductivity evaluator").get<std::string>("height key") != getKey(domain_,"ponded_depth"));

  // Now continue as usual for overland head
  OverlandPressureFlow::SetupPhysicalEvaluators_(S);
}

} // namespace
} // namespace
