/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for overland flow.
------------------------------------------------------------------------- */

#include "energy_surface_ice.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<EnergySurfaceIce> EnergySurfaceIce::reg_("surface energy");

} // namespace
} // namespace
#include "surface_ice_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<FieldEvaluator,SurfaceIceEnergyEvaluator> SurfaceIceEnergyEvaluator::reg_("surface ice energy");

} //namespace
} //namespace
#include "surface_subgrid_ice_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<FieldEvaluator,SurfaceSubgridIceEnergyEvaluator> SurfaceSubgridIceEnergyEvaluator::reg_("surface subgrid ice energy");

} //namespace
} //namespace
