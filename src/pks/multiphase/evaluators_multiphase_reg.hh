/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Process Kernels

  Miscalleneous general purpose evaluators for various flow models.
*/

#include "MassDensityGas.hh"
#include "MoleFractionLiquid.hh"
#include "NCP_F.hh"
#include "NCP_HenryLaw.hh"
#include "NCP_MolarDensities.hh"
#include "NCP_MoleFractions.hh"
#include "ProductEvaluator.hh"
#include "SaturationEvaluator.hh"
#include "TotalComponentStorage.hh"

namespace Amanzi {
namespace Multiphase {

Utils::RegisteredFactory<Evaluator, MassDensityGas> MassDensityGas::fac_("mass density gas");
Utils::RegisteredFactory<Evaluator, MoleFractionLiquid> MoleFractionLiquid::fac_("mole fraction liquid");

Utils::RegisteredFactory<Evaluator, NCP_F> NCP_F::fac_("ncp saturation");
Utils::RegisteredFactory<Evaluator, NCP_HenryLaw> NCP_HenryLaw::fac_("ncp henry law");
Utils::RegisteredFactory<Evaluator, NCP_MolarDensities> NCP_MolarDensities::fac_("ncp molar densities");
Utils::RegisteredFactory<Evaluator, NCP_MoleFractions> NCP_MoleFractions::fac_("ncp mole fraction");

Utils::RegisteredFactory<Evaluator, ProductEvaluator> ProductEvaluator::fac_("product");
Utils::RegisteredFactory<Evaluator, SaturationEvaluator> SaturationEvaluator::fac_("saturation");

Utils::RegisteredFactory<Evaluator, TotalComponentStorage> TotalComponentStorage::fac_("storage component");

} // namespace Multiphase
} // namespace Amanzi
