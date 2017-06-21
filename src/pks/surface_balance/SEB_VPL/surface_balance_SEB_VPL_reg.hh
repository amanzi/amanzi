/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon, Adam Atchley, Satish Karra
 *
 * ------------------------------------------------------------------------- */


#include "surface_balance_evaluator_VPL.hh"
#include "surface_balance_SEB_VPL.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory<SurfaceBalanceSEBVPL>
SurfaceBalanceSEBVPL::reg_("surface balance SEB VPL");

Utils::RegisteredFactory<FieldEvaluator,SurfaceBalanceEvaluatorVPL>
SurfaceBalanceEvaluatorVPL::reg_("surface balance SEB VPL");

} // namespace
} // namespace
