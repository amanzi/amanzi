/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Self-registering factory of EOS models.
*/

#include "EOSDensityEvaluator.hh"
#include "EOSDensityFactory.hh"
#include "EOS_DensityIdealGas.hh"
#include "EOS_DensityVaporInGas.hh"
#include "EOS_DensityWater.hh"
#include "EOS_DensityWaterFEHM.hh"
#include "EOS_DensityWaterTabular.hh"
#include "IsobaricEOSEvaluator.hh"
#include "MolarFractionGasEvaluator.hh"
#include "SaturatedVaporPressure_Water.hh"
#include "SaturatedVaporPressureFactory.hh"
#include "ViscosityConstant.hh"
#include "ViscosityWater.hh"
#include "ViscosityWaterFEHM.hh"
#include "ViscosityWaterTabular.hh"
#include "ViscosityBaseFactory.hh"
#include "ViscosityEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<FieldEvaluator, EOSDensityEvaluator> EOSDensityEvaluator::factory_("eos");
Utils::RegisteredFactory<FieldEvaluator, IsobaricEOSEvaluator> IsobaricEOSEvaluator::factory_("isobaric eos");
Utils::RegisteredFactory<FieldEvaluator, MolarFractionGasEvaluator> MolarFractionGasEvaluator::factory_("molar fraction gas");
Utils::RegisteredFactory<FieldEvaluator, ViscosityEvaluator> ViscosityEvaluator::factory_("viscosity");

Utils::RegisteredFactory<EOS_Density, EOS_DensityIdealGas> EOS_DensityIdealGas::factory_("ideal gas");
Utils::RegisteredFactory<EOS_Density, EOS_DensityVaporInGas> EOS_DensityVaporInGas::factory_("vapor in gas");
Utils::RegisteredFactory<EOS_Density, EOS_DensityWater> EOS_DensityWater::factory_("liquid water 0-30C");
Utils::RegisteredFactory<EOS_Density, EOS_DensityWaterFEHM> EOS_DensityWaterFEHM::factory_("liquid water FEHM");
Utils::RegisteredFactory<EOS_Density, EOS_DensityWaterTabular> EOS_DensityWaterTabular::factory_("liquid water tabular");

Utils::RegisteredFactory<SaturatedVaporPressure, SaturatedVaporPressure_Water> SaturatedVaporPressure_Water::factory_("water vapor over water/ice");

Utils::RegisteredFactory<ViscosityBase, ViscosityConstant> ViscosityConstant::factory_("constant");
Utils::RegisteredFactory<ViscosityBase, ViscosityWater> ViscosityWater::factory_("liquid water 0-30C");
Utils::RegisteredFactory<ViscosityBase, ViscosityWaterFEHM> ViscosityWaterFEHM::factory_("liquid water FEHM");
Utils::RegisteredFactory<ViscosityBase, ViscosityWaterTabular> ViscosityWaterTabular::factory_("liquid water tabular");

}  // namespace AmanziEOS
}  // namespace Amanzi

