/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  A field evaluator for an unchanging cell volume.
*/

#include "EvaluatorIndependentFromFile.hh"
#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorIndependentConstant.hh"
#include "EvaluatorIndependentTensorFunction.hh"
#include "EvaluatorIndependentPatchFunction.hh"
#include "EvaluatorSecondaryMonotypeFromFunction.hh"
#include "EvaluatorSecondaryMonotypeAdditive.hh"
#include "EvaluatorSecondaryMonotypeMultiplicative.hh"
#include "EvaluatorPrimaryStaticMesh.hh"
#include "EvaluatorSecondaryMeshedQuantity.hh"
#include "EvaluatorSecondaryVectorAsPatch.hh"
#include "EvaluatorAggregateBCs.hh"

#include "registration_macro.hh"

namespace Amanzi {

REGISTER(EvaluatorIndependentFunction);
REGISTER(EvaluatorIndependentFromFile);
REGISTER(EvaluatorIndependentConstant);
REGISTER(EvaluatorIndependentTensorFunction);
REGISTER(EvaluatorIndependentPatchFunction);
REGISTER(EvaluatorSecondaryMonotypeFromFunction);
template <>
REGISTER(EvaluatorSecondaryMonotypeAdditiveCV);
template <>
REGISTER(EvaluatorSecondaryMonotypeMultiplicativeCV);

const std::string EvaluatorPrimaryStaticMesh::eval_type = "static mesh";
REGISTER(EvaluatorPrimaryStaticMesh);

template <>
REGISTER(EvaluatorCellVolume);
template <>
REGISTER(EvaluatorMeshElevation);
template <>
REGISTER(EvaluatorMeshSlopeMagnitude);
template <>
REGISTER(EvaluatorMeshAspect);

REGISTER(EvaluatorAggregateBCs);
REGISTER(EvaluatorSecondaryVectorAsPatch);

} // namespace Amanzi
