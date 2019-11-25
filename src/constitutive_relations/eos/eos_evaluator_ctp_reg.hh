/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#include "eos_evaluator_ctp.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EOSEvaluatorCTP> EOSEvaluatorCTP::factory_("eos ctp");

}
}
