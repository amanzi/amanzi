/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

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
#include "H2O_ThermalConductivity.hh"
#include "H2O_Viscosity.hh"
#include "H2O_ViscosityFEHM.hh"
#include "H2O_ViscosityTabular.hh"
#include "IdealGas_Density.hh"
#include "IsobaricEOSEvaluator.hh"
#include "MolarFractionGasEvaluator.hh"
#include "ThermalConductivityConstant.hh"
#include "VaporInGas_Density.hh"
#include "VaporInGas_Diffusion.hh"
#include "ViscosityConstant.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<Evaluator, EOSDensityEvaluator> EOSDensityEvaluator::reg_("eos");
Utils::RegisteredFactory<Evaluator, IsobaricEOSEvaluator>
  IsobaricEOSEvaluator::reg_("isobaric eos");
Utils::RegisteredFactory<Evaluator, MolarFractionGasEvaluator>
  MolarFractionGasEvaluator::reg_("molar fraction gas");
Utils::RegisteredFactory<Evaluator, EOSViscosityEvaluator> EOSViscosityEvaluator::reg_("viscosity");

Utils::RegisteredFactory<EOS_Density, IdealGas_Density> IdealGas_Density::reg_("ideal gas");
Utils::RegisteredFactory<EOS_Density, VaporInGas_Density> VaporInGas_Density::reg_("vapor in gas");
Utils::RegisteredFactory<EOS_Density, H2O_Density> H2O_Density::reg_("liquid water 0-30C");
Utils::RegisteredFactory<EOS_Density, H2O_DensityFEHM> H2O_DensityFEHM::reg_("liquid water FEHM");
Utils::RegisteredFactory<EOS_Density, H2O_DensityTabular>
  H2O_DensityTabular::reg_("liquid water tabular");

Utils::RegisteredFactory<EOS_SaturatedVaporPressure, H2O_SaturatedVaporPressure>
  H2O_SaturatedVaporPressure::reg_("water vapor over water/ice");

Utils::RegisteredFactory<EOS_Viscosity, H2O_Viscosity> H2O_Viscosity::reg_("liquid water 0-30C");
Utils::RegisteredFactory<EOS_Viscosity, ViscosityConstant> ViscosityConstant::reg_("constant");
Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityFEHM>
  H2O_ViscosityFEHM::reg_("liquid water FEHM");
Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityTabular>
  H2O_ViscosityTabular::reg_("liquid water tabular");

Utils::RegisteredFactory<EOS_Diffusion, VaporInGas_Diffusion>
  VaporInGas_Diffusion::reg_("vapor in gas");

Utils::RegisteredFactory<COM_Tortuosity, COM_MillingtonQuirk>
  COM_MillingtonQuirk::reg_("Millington Quirk");

Utils::RegisteredFactory<EOS_ThermalConductivity, H2O_ThermalConductivity>
  H2O_ThermalConductivity::reg_("liquid water");
Utils::RegisteredFactory<EOS_ThermalConductivity, ThermalConductivityConstant>
  ThermalConductivityConstant::reg_("constant");

} // namespace AmanziEOS
} // namespace Amanzi
