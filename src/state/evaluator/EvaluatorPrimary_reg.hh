/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#include "evaluator/EvaluatorPrimary.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator,EvaluatorPrimary> EvaluatorPrimary::fac_("primary variable");

} // namespace
