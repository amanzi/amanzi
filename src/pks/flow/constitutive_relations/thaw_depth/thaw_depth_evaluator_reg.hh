/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "thaw_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ThawDepthEvaluator> ThawDepthEvaluator::reg_("thaw depth");

}
}
