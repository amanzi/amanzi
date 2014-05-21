/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "iem_quadratic.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

Utils::RegisteredFactory<IEM,IEMQuadratic> IEMQuadratic::factory_("quadratic");

} // namespace
} // namespace
} // namespace
