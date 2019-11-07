/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid).
*/

#include "EffectivePressureEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

// registry of method
Utils::RegisteredFactory<FieldEvaluator, EffectivePressureEvaluator>
    EffectivePressureEvaluator::factory_("effective_pressure");

}  // namespace AmanziEOS
}  // namespace Amanzi


