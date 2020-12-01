#include "water_table_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WaterTableDepthEvaluator> WaterTableDepthEvaluator::reg_("water table depth");

}
}
