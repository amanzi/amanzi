/*
  Process Kernels
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscalleneous general purpose evaluators for various flow models.
*/

#include "MoleFractionLiquid.hh"
#include "NCP_HenryLaw.hh"
#include "NCP_MolarDensities.hh"
#include "NCP_MoleFractions.hh"
#include "ProductEvaluator.hh"
#include "SaturationGasEvaluator.hh"
#include "TccGas.hh"
#include "TccLiquid.hh"
#include "TotalComponentStorage.hh"
#include "TotalComponentStorage_Tcc.hh"
#include "TotalComponentStorage_Jaffre.hh"
#include "TotalWaterStorage.hh"

namespace Amanzi {
namespace Multiphase {

Utils::RegisteredFactory<Evaluator, MoleFractionLiquid>
  MoleFractionLiquid::fac_("mole fraction liquid");
Utils::RegisteredFactory<Evaluator, NCP_HenryLaw> NCP_HenryLaw::fac_("ncp henry law");
Utils::RegisteredFactory<Evaluator, NCP_MolarDensities>
  NCP_MolarDensities::fac_("ncp molar densities");
Utils::RegisteredFactory<Evaluator, NCP_MoleFractions>
  NCP_MoleFractions::fac_("ncp mole fraction gas");
Utils::RegisteredFactory<Evaluator, ProductEvaluator> ProductEvaluator::fac_("product");
Utils::RegisteredFactory<Evaluator, SaturationGasEvaluator>
  SaturationGasEvaluator::fac_("saturation gas");
Utils::RegisteredFactory<Evaluator, TotalComponentStorage>
  TotalComponentStorage::fac_("storage component");
Utils::RegisteredFactory<Evaluator, TotalComponentStorage_Tcc>
  TotalComponentStorage_Tcc::fac_("storage component tcc");
Utils::RegisteredFactory<Evaluator, TotalComponentStorage_Jaffre>
  TotalComponentStorage_Jaffre::fac_("storage component jaffre");
Utils::RegisteredFactory<Evaluator, TotalWaterStorage> TotalWaterStorage::fac_("storage water");
Utils::RegisteredFactory<Evaluator, TccGas> TccGas::fac_("tcc gas");
Utils::RegisteredFactory<Evaluator, TccLiquid> TccLiquid::fac_("tcc liquid");

} // namespace Multiphase
} // namespace Amanzi
