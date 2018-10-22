#include "seb_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,SEBEvaluator> SEBEvaluator::reg_("surface balance");

} //namespace
} //namespace
