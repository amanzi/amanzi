/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  MPC PK

  Registration of MPC evalautors.
*/

#include "HeatDiffusionMatrixFracture.hh"
#include "SoluteDiffusionMatrixFracture.hh"
#include "WaterStorageStressSplit.hh"

namespace Amanzi {

Utils::RegisteredFactory<Evaluator, SoluteDiffusionMatrixFracture>
  SoluteDiffusionMatrixFracture::reg_("solute diffusion to matrix");

Utils::RegisteredFactory<Evaluator, HeatDiffusionMatrixFracture> HeatDiffusionMatrixFracture::reg_(
  "heat diffusion to matrix");

Utils::RegisteredFactory<Evaluator, WaterStorageStressSplit> WaterStorageStressSplit::reg_(
  "water storage stress split");

} // namespace Amanzi
