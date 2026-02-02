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
#include "ThermodynamicStateEvaluators.hh"
#include "VolumetricStrainEvaluator.hh"

namespace Amanzi {
namespace Evaluators {

Utils::RegisteredFactory<Evaluator, PorosityEvaluator> PorosityEvaluator::reg_("porosity");
Utils::RegisteredFactory<Evaluator, VolumetricStrainEvaluator> VolumetricStrainEvaluator::reg_(
  "volumetric strain");

Utils::RegisteredFactory<Evaluator, ThermodynamicStateEvaluator> ThermodynamicStateEvaluator::reg_("iapws97 state");
Utils::RegisteredFactory<Evaluator, DensityEvaluator> DensityEvaluator::reg_("iapws97 density");
Utils::RegisteredFactory<Evaluator, TemperatureEvaluator> TemperatureEvaluator::reg_("iapws97 temperature");
Utils::RegisteredFactory<Evaluator, ViscosityEvaluator> ViscosityEvaluator::reg_("iapws97 viscosity");
Utils::RegisteredFactory<Evaluator, IsothermalCompressibilityEvaluator>
   IsothermalCompressibilityEvaluator::reg_("iapws97 isothermal compressibility");

} // namespace Evaluators
} // namespace Amanzi
