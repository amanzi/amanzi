#include "evaluator/EvaluatorSecondaryMeshedQuantity.hh"

namespace Amanzi {

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorCellVolume> EvaluatorCellVolume::fac_("cell volume");
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorFaceArea> EvaluatorFaceArea::fac_("face area");
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorMeshElevation> EvaluatorMeshElevation::fac_("meshed elevation");
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorMeshSlopeMagnitude> EvaluatorMeshSlopeMagnitude::fac_("meshed slope magnitude");
}
