/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "elevation_evaluator_column.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ElevationEvaluatorColumn> ElevationEvaluatorColumn::reg_("elevation evaluator");

}
}
