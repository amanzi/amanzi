/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Ethan Coon, Adam Atchley, Satish Karra
 *
 * ------------------------------------------------------------------------- */


#include "surface_balance_evaluator.hh"
#include "surface_balance_SEB.hh"

namespace Amanzi {
namespace SurfaceBalance {

RegisteredPKFactory_ATS<SurfaceBalanceSEB>
SurfaceBalanceSEB::reg_("surface balance SEB");

Utils::RegisteredFactory<FieldEvaluator,SurfaceBalanceEvaluator>
SurfaceBalanceEvaluator::reg_("surface balance SEB");

} // namespace
} // namespace
