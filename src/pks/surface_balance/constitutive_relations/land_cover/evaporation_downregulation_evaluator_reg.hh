#include "evaporation_downregulation_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<FieldEvaluator,EvaporationDownregulationEvaluator> EvaporationDownregulationEvaluator::reg_("evaporation downregulation via soil resistance");

} //namespace
} //namespace
} //namespace
