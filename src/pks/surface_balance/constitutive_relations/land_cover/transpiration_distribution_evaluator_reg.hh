#include "transpiration_distribution_evaluator.hh"

namespace Amanzi {
namespace LandCover {
namespace Relations {

Utils::RegisteredFactory<FieldEvaluator,TranspirationDistributionEvaluator> TranspirationDistributionEvaluator::reg_("transpiration distribution via rooting depth");

} //namespace
} //namespace
} //namespace
