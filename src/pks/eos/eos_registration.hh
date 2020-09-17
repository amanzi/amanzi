/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Self-registering factory of EOS models.
*/

#include "EOSFactory.hh"
#include "EOSEvaluator.hh"
#include "EOS_Constant.hh"
#include "EOS_IdealGas.hh"
#include "EOS_VaporInGas.hh"
#include "EOS_Water.hh"
#include "IsobaricEOSEvaluator.hh"
#include "MolarFractionGasEvaluator.hh"
#include "SaturatedVaporPressure_Water.hh"
#include "SaturatedVaporPressureFactory.hh"
#include "Viscosity_Constant.hh"
#include "Viscosity_Water.hh"
#include "ViscosityBaseFactory.hh"
#include "ViscosityEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<FieldEvaluator, EOSEvaluator> EOSEvaluator::factory_("eos");
Utils::RegisteredFactory<FieldEvaluator, IsobaricEOSEvaluator> IsobaricEOSEvaluator::factory_("isobaric eos");
Utils::RegisteredFactory<FieldEvaluator, MolarFractionGasEvaluator> MolarFractionGasEvaluator::factory_("molar fraction gas");
Utils::RegisteredFactory<FieldEvaluator, ViscosityEvaluator> ViscosityEvaluator::factory_("viscosity");

Utils::RegisteredFactory<EOS, EOS_Constant> EOS_Constant::factory_("constant");
Utils::RegisteredFactory<EOS, EOS_IdealGas> EOS_IdealGas::factory_("ideal gas");
Utils::RegisteredFactory<EOS, EOS_VaporInGas> EOS_VaporInGas::factory_("vapor in gas");
Utils::RegisteredFactory<EOS, EOS_Water> EOS_Water::factory_("liquid water");

Utils::RegisteredFactory<SaturatedVaporPressure, SaturatedVaporPressure_Water> SaturatedVaporPressure_Water::factory_("water vapor over water/ice");

Utils::RegisteredFactory<Viscosity_Base, Viscosity_Constant> Viscosity_Constant::factory_("constant");
Utils::RegisteredFactory<Viscosity_Base, Viscosity_Water> Viscosity_Water::factory_("liquid water");

}  // namespace AmanziEOS
}  // namespace Amanzi

