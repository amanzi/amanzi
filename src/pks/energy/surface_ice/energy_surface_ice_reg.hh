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

RegisteredPKFactory_ATS<EnergySurfaceIce> EnergySurfaceIce::reg_("surface energy");

} // namespace
} // namespace
