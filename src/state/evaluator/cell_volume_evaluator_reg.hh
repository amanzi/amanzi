/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging cell volume.

------------------------------------------------------------------------- */

#include "cell_volume_evaluator.hh"

namespace Amanzi {

Utils::RegisteredFactory<Evaluator, CellVolumeEvaluator>
    CellVolumeEvaluator::fac_("cell volume");

} // namespace Amanzi
