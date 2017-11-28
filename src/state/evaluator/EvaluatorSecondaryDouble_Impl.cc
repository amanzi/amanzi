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

#include "EvaluatorSecondary.hh"


namespace Amanzi {

// ---------------------------------------------------------------------------
// Updates the derivative for doubles
// ---------------------------------------------------------------------------
template<>
void EvaluatorSecondary<double>::UpdateDerivative_(State& S, const Key& wrt_key) {
  Key dmy_key = std::string("d")+my_key_+std::string("_d")+wrt_key;
  if (!S.HasData(dmy_key, my_tag_)) {
    // or create the field.  Note we have to do extra work that is normally
    // done by State in setup.
    S.Require<double>(dmy_key, my_tag_, my_key_);
    auto dmy_ptr = Teuchos::rcp(new double);
    *dmy_ptr = 0.;
    S.SetPtr(dmy_key, my_tag_, my_key_, dmy_ptr);
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_initialized();
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_io_vis(false);
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_io_checkpoint(false);
  }

  double& dmy = S.GetW<double>(dmy_key, my_tag_, my_key_);
  dmy = 0.;

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {

    
    if (wrt_key == *dep) {
      // partial F / partial x
      double tmp;
      EvaluatePartialDerivative_(S, wrt_key, tmp);
      dmy += tmp;
      // std::cout << "  chain rule: p" << my_key_ << "_p" << wrt_key << " = " << *tmp << std::endl;

    } else if (S.GetEvaluator(*dep)->IsDependency(S, wrt_key)) {
      // partial F / partial dep * ddep/dx
      // -- ddep/dx
      Key ddep_key = std::string("d")+*dep+std::string("_d")+wrt_key;
      auto& ddep = S.Get<double>(ddep_key);

      // -- partial F / partial dep
      double tmp;
      EvaluatePartialDerivative_(S, *dep, tmp);

      // std::cout << "  chain rule: p" << my_key_ << "_p" << *dep << " * p" << *dep << "_p" << wrt_key << " = "
      // 	<< ddep << " * " << *tmp << std::endl;
      dmy += ddep * tmp;
    }
  }
}

} // namespace
