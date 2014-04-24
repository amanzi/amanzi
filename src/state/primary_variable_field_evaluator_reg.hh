/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,PrimaryVariableFieldEvaluator> PrimaryVariableFieldEvaluator::fac_("primary variable");

} // namespace
