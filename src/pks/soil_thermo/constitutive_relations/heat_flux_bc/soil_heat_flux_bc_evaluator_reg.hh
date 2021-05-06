#include "soil_heat_flux_bc_evaluator.hh"

namespace Amanzi {
namespace SoilThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SoilHeatFluxBCEvaluator> SoilHeatFluxBCEvaluator::factory_("heat flux bc");

} //namespace
} //namespace
