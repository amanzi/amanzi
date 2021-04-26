/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "carbon_decomposition_rate_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,CarbonDecomposeRateEvaluator> CarbonDecomposeRateEvaluator::reg_("carbon decomposition rate");

}
}
