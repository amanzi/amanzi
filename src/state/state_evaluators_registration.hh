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

#include "EvaluatorCellVolume.hh"
#include "EvaluatorDeformingCellVolume.hh"
#include "EvaluatorIndependentFromFile.hh"
#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorIndependentConstant.hh"
#include "EvaluatorIndependentTensorFunction.hh"
#include "EvaluatorMultiplicativeReciprocal.hh"
#include "EvaluatorSecondaryMonotypeFromFunction.hh"
#include "EvaluatorPrimary.hh"
#include "EvaluatorTemporalInterpolation.hh"
#include "EvaluatorVelocityReconstruction.hh"

namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorCellVolume> EvaluatorCellVolume::fac_("cell volume");
Utils::RegisteredFactory<Evaluator, EvaluatorDeformingCellVolume>
  EvaluatorDeformingCellVolume::fac_("deforming cell volume");

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction>
  EvaluatorIndependentFunction::fac_("independent variable");
Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFromFile>
  EvaluatorIndependentFromFile::fac_("independent variable from file");
Utils::RegisteredFactory<Evaluator, EvaluatorIndependentConstant>
  EvaluatorIndependentConstant::fac_("independent variable constant");

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFunction>
  EvaluatorIndependentTensorFunction::fac_("independent variable tensor");

Utils::RegisteredFactory<Evaluator, EvaluatorMultiplicativeReciprocal>
  EvaluatorMultiplicativeReciprocal::fac_("multiplicative reciprocal");

Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeFromFunction>
  EvaluatorSecondaryMonotypeFromFunction::fac_("secondary variable from function");

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorPrimaryCV> EvaluatorPrimaryCV::fac_(
  "primary variable");

Utils::RegisteredFactory<Evaluator, EvaluatorTemporalInterpolation>
  EvaluatorTemporalInterpolation::fac_("temporal interpolation");

Utils::RegisteredFactory<Evaluator, EvaluatorVelocityReconstruction>
  EvaluatorVelocityReconstruction::fac_("velocity reconstruction");

} // namespace Amanzi
