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
template<>
Utils::RegisteredFactory<Evaluator,EvaluatorPrimary<double>>
EvaluatorPrimary<double>::fac_("primary variable double");

template<>
Utils::RegisteredFactory<Evaluator,EvaluatorPrimary<CompositeVector,CompositeVectorSpace>>
EvaluatorPrimary<CompositeVector,CompositeVectorSpace>::fac_("primary variable");

} // namespace
