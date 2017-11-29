/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator base for secondary variables.  A Evaluator is a node in the
Phalanx-like dependency tree.

Secondary variable evaluators, such as equations of state, water retention evaluators,
internal energy evaluators, etc should inherit this class, implementing the
missing Update_() and UpdateFieldDerivative_() methods.

Algebraic secondary evaluator where all dependencies and this are
doubles.

------------------------------------------------------------------------- */

template<>
void EvaluatorSecondary<double>::UpdateDerivative_(State& S, const Key& wrt_key);
