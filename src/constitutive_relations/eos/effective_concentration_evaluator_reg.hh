/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EffectiveConcentrationEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "effective_concentration_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EffectiveConcentrationEvaluator> EffectiveConcentrationEvaluator::factory_("effective_concentration");

}
}


