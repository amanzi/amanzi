/*
  Process Kernels
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Miscalleneous general purpose evaluators for various flow models.
*/

#include "SaturationGasEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

Utils::RegisteredFactory<Evaluator, SaturationGasEvaluator> SaturationGasEvaluator::fac_("saturation gas");

}  // namespace Multiphase
}  // namespace Amanzi

