/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Evaluators

  Self-registering factory for common evaluators.
*/

#include "PorosityEvaluator.hh"
#include "VolumetricStrainEvaluator.hh"

namespace Amanzi {
namespace Evaluators {

Utils::RegisteredFactory<Evaluator, PorosityEvaluator> PorosityEvaluator::reg_("porosity");
Utils::RegisteredFactory<Evaluator, VolumetricStrainEvaluator> VolumetricStrainEvaluator::reg_(
  "volumetric strain");

} // namespace Evaluators
} // namespace Amanzi
