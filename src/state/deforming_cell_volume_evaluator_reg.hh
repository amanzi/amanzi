/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an changing cell volume.

------------------------------------------------------------------------- */

#include "deforming_cell_volume_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<FieldEvaluator,DeformingCellVolumeEvaluator>
DeformingCellVolumeEvaluator::fac_("deforming cell volume");

} // namespace
