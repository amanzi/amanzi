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

#include "IAPWS95_StateEvaluators.hh"
#include "PorosityEvaluator.hh"
#include "IAPWS97_StateEvaluators.hh"
#include "VolumetricStrainEvaluator.hh"

namespace Amanzi {
namespace Evaluators {

Utils::RegisteredFactory<Evaluator, PorosityEvaluator> PorosityEvaluator::reg_("porosity");
Utils::RegisteredFactory<Evaluator, VolumetricStrainEvaluator> VolumetricStrainEvaluator::reg_(
  "volumetric strain");

Utils::RegisteredFactory<Evaluator, IAPWS95_StateEvaluator> IAPWS95_StateEvaluator::reg_("iapws95 state");
Utils::RegisteredFactory<Evaluator, IAPWS95_DensityEvaluator> IAPWS95_DensityEvaluator::reg_("iapws95 density");
Utils::RegisteredFactory<Evaluator, IAPWS95_ViscosityEvaluator> IAPWS95_ViscosityEvaluator::reg_("iapws95 viscosity");
Utils::RegisteredFactory<Evaluator, IAPWS95_InternalEnergyEvaluator> IAPWS95_InternalEnergyEvaluator::reg_("iapws95 internal energy");

Utils::RegisteredFactory<Evaluator, IAPWS97_StateEvaluator> IAPWS97_StateEvaluator::reg_("iapws97 state");
Utils::RegisteredFactory<Evaluator, IAPWS97_DensityEvaluator> IAPWS97_DensityEvaluator::reg_("iapws97 density");
Utils::RegisteredFactory<Evaluator, IAPWS97_TemperatureEvaluator> IAPWS97_TemperatureEvaluator::reg_("iapws97 temperature");
Utils::RegisteredFactory<Evaluator, IAPWS97_ViscosityEvaluator> IAPWS97_ViscosityEvaluator::reg_("iapws97 viscosity");
Utils::RegisteredFactory<Evaluator, IAPWS97_IsothermalCompressibilityEvaluator>
   IAPWS97_IsothermalCompressibilityEvaluator::reg_("iapws97 isothermal compressibility");

} // namespace Evaluators
} // namespace Amanzi
