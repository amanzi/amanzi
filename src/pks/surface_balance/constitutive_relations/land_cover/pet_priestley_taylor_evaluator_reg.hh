#include "pet_preistley_taylor_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<FieldEvaluator,PETPriestleyTaylorEvaluator> PETPriestleyTaylorEvaluator::reg_("potential evapotranspiration, Priestley-Taylor");

} //namespace
} //namespace
} //namespace
