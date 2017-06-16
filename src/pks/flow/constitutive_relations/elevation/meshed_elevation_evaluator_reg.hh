/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "meshed_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,MeshedElevationEvaluator> MeshedElevationEvaluator::reg_("meshed elevation");

}
}
