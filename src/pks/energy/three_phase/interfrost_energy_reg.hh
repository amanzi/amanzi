/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "interfrost_energy.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<InterfrostEnergy> InterfrostEnergy::reg_("interfrost energy");

} // namespace Energy
} // namespace Amanzi
