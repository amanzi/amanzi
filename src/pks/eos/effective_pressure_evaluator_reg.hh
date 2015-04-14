/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.
*/

#include "effective_pressure_evaluator.hh"

namespace Amanzi {
namespace EOS {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EffectivePressureEvaluator> EffectivePressureEvaluator::factory_("effective_pressure");

}  // namespace EOS
}  // namespace Amanzi


