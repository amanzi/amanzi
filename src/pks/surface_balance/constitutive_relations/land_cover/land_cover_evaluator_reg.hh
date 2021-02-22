/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#include "drainage_evaluator.hh"
#include "interception_fraction_evaluator.hh"
#include "pet_priestley_taylor_evaluator.hh"
#include "plant_wilting_factor_evaluator.hh"
#include "rooting_depth_fraction_evaluator.hh"
#include "snow_meltrate_evaluator.hh"
#include "transpiration_distribution_evaluator.hh"
#include "radiation_balance_evaluator.hh"
#include "albedo_evaluator.hh"
#include "area_fractions_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,DrainageEvaluator> DrainageEvaluator::factory_("canopy drainage");

Utils::RegisteredFactory<FieldEvaluator,InterceptionFractionEvaluator> InterceptionFractionEvaluator::reg_("interception fraction");

Utils::RegisteredFactory<FieldEvaluator,PETPriestleyTaylorEvaluator> PETPriestleyTaylorEvaluator::reg_("potential evapotranspiration, Priestley-Taylor");

Utils::RegisteredFactory<FieldEvaluator,PlantWiltingFactorEvaluator> PlantWiltingFactorEvaluator::reg_("plant wilting factor");

Utils::RegisteredFactory<FieldEvaluator,RootingDepthFractionEvaluator> RootingDepthFractionEvaluator::reg_("rooting depth fraction");

Utils::RegisteredFactory<FieldEvaluator,SnowMeltRateEvaluator> SnowMeltRateEvaluator::reg_("snow melt rate");

Utils::RegisteredFactory<FieldEvaluator,TranspirationDistributionEvaluator> TranspirationDistributionEvaluator::reg_("transpiration distribution via rooting depth");

Utils::RegisteredFactory<FieldEvaluator,RadiationBalanceEvaluator> RadiationBalanceEvaluator::reg_("radiation balance, surface and canopy");

Utils::RegisteredFactory<FieldEvaluator,AlbedoEvaluator> AlbedoEvaluator::reg_("ground albedo");

Utils::RegisteredFactory<FieldEvaluator,AreaFractionsEvaluator> AreaFractionsEvaluator::reg_("snow area fraction");


} // namespace
} // namespace
} // namespace
