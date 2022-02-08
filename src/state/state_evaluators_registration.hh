/*
  State

  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator for an unchanging cell volume.
*/

#include "EvaluatorCellVolume.hh"
#include "EvaluatorDeformingCellVolume.hh"
#include "EvaluatorIndependentFromFile.hh"
#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorIndependentConstant.hh"
#include "EvaluatorMultiplicativeReciprocal.hh"
#include "EvaluatorSecondaryMonotypeFromFunction.hh"

namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorCellVolume> EvaluatorCellVolume::fac_("cell volume");
Utils::RegisteredFactory<Evaluator, EvaluatorDeformingCellVolume> EvaluatorDeformingCellVolume::fac_("deforming cell volume");

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction> EvaluatorIndependentFunction::fac_("independent variable");
Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFromFile> EvaluatorIndependentFromFile::fac_("independent variable from file");
Utils::RegisteredFactory<Evaluator, EvaluatorIndependentConstant> EvaluatorIndependentConstant::fac_("independent variable constant");

Utils::RegisteredFactory<Evaluator, EvaluatorMultiplicativeReciprocal> EvaluatorMultiplicativeReciprocal::fac_("multiplicative reciprocal");

Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeFromFunction> EvaluatorSecondaryMonotypeFromFunction::fac_("secondary variable from function");

} // namespace

