/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<AdvectionDiffusion> AdvectionDiffusion::reg_("advection-diffusion energy");

} // namespace Energy
} // namespace Amanzi
