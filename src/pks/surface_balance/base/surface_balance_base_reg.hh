/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon, Adam Atchley, Satish Karra
 *
 * ------------------------------------------------------------------------- */


#include "surface_balance_base.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory_ATS<SurfaceBalanceBase>
SurfaceBalanceBase::reg_("general surface balance");

} // namespace
} // namespace
