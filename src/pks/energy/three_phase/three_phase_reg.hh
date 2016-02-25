/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "three_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory_ATS<ThreePhase> ThreePhase::reg_("three-phase energy");

} // namespace Energy
} // namespace Amanzi
