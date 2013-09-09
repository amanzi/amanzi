/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging cell volume.

------------------------------------------------------------------------- */

#include "cell_volume_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,CellVolumeEvaluator> 
CellVolumeEvaluator::fac_("cell volume");

}
