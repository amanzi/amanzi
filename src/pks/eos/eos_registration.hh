/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Self-registering factory of EOS models.
*/

#include "COM_MillingtonQuirk.hh"

#include "EOSDensityEvaluator.hh"
#include "EOSFactory.hh"
#include "EOSViscosityEvaluator.hh"
#include "H2O_DensityTabular.hh"
#include "H2O_Density.hh"
#include "H2O_DensityFEHM.hh"
#include "H2O_SaturatedVaporPressure.hh"
#include "H2O_Viscosity.hh"
#include "H2O_ViscosityConstant.hh"
#include "H2O_ViscosityFEHM.hh"
#include "H2O_ViscosityTabular.hh"
#include "IdealGas_Density.hh"
#include "IsobaricEOSEvaluator.hh"
#include "MolarFractionGasEvaluator.hh"
#include "VaporInGas_Density.hh"
#include "VaporInGas_Diffusion.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<Evaluator, EOSDensityEvaluator> EOSDensityEvaluator::factory_("eos");
Utils::RegisteredFactory<Evaluator, IsobaricEOSEvaluator>
  IsobaricEOSEvaluator::factory_("isobaric eos");
Utils::RegisteredFactory<Evaluator, MolarFractionGasEvaluator>
  MolarFractionGasEvaluator::factory_("molar fraction gas");
Utils::RegisteredFactory<Evaluator, EOSViscosityEvaluator>
  EOSViscosityEvaluator::factory_("viscosity");

Utils::RegisteredFactory<EOS_Density, IdealGas_Density> IdealGas_Density::factory_("ideal gas");
Utils::RegisteredFactory<EOS_Density, VaporInGas_Density>
  VaporInGas_Density::factory_("vapor in gas");
Utils::RegisteredFactory<EOS_Density, H2O_Density> H2O_Density::factory_("liquid water 0-30C");
Utils::RegisteredFactory<EOS_Density, H2O_DensityFEHM>
  H2O_DensityFEHM::factory_("liquid water FEHM");
Utils::RegisteredFactory<EOS_Density, H2O_DensityTabular>
  H2O_DensityTabular::factory_("liquid water tabular");

Utils::RegisteredFactory<EOS_SaturatedVaporPressure, H2O_SaturatedVaporPressure>
  H2O_SaturatedVaporPressure::factory_("water vapor over water/ice");

Utils::RegisteredFactory<EOS_Viscosity, H2O_Viscosity>
  H2O_Viscosity::factory_("liquid water 0-30C");
Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityConstant>
  H2O_ViscosityConstant::factory_("constant");
Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityFEHM>
  H2O_ViscosityFEHM::factory_("liquid water FEHM");
Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityTabular>
  H2O_ViscosityTabular::factory_("liquid water tabular");

Utils::RegisteredFactory<EOS_Diffusion, VaporInGas_Diffusion>
  VaporInGas_Diffusion::factory_("vapor in gas");

Utils::RegisteredFactory<COM_Tortuosity, COM_MillingtonQuirk>
  COM_MillingtonQuirk::factory_("Millington Quirk");

} // namespace AmanziEOS
} // namespace Amanzi
