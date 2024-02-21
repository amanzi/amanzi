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

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentPatchFunction>
  EvaluatorIndependentPatchFunction::fac_("independent variable patch");

Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeFromFunction>
  EvaluatorSecondaryMonotypeFromFunction::fac_("secondary variable from function");

template <>
REGISTER(EvaluatorSecondaryMonotypeAdditiveCV);

Utils::RegisteredFactory<Evaluator, EvaluatorPrimaryStaticMesh>
  EvaluatorPrimaryStaticMesh::fac_("static mesh");


template <>
REGISTER(EvaluatorCellVolume);
template <>
REGISTER(EvaluatorMeshElevation);
template <>
REGISTER(EvaluatorMeshSlopeMagnitude);

Utils::RegisteredFactory<Evaluator, EvaluatorAggregateBCs>
  EvaluatorAggregateBCs::fac_("boundary condition aggregator");



} // namespace Amanzi
