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
CompositeVectors.

------------------------------------------------------------------------- */


// ---------------------------------------------------------------------------
// Ensures that dependencies provide the vector structure we need for this.
// ---------------------------------------------------------------------------
template<>
void
EvaluatorSecondary<CompositeVector, CompositeVectorSpace>::EnsureCompatibility(State& S) {
  // Ensure my field exists.  Requirements should be already set.
  ASSERT(my_key_ != std::string(""));
  CompositeVectorSpace dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key_, my_tag_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S.GetRecordW(my_key_, my_tag_, my_key_).set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (dep_fac.Mesh() != Teuchos::null) {
    dep_fac.SetOwned(false);

    // Loop over my dependencies, ensuring they meet the requirements.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      auto& fac = S.Require<CompositeVector,CompositeVectorSpace>(*key, my_tag_);
      fac.Update(dep_fac);
    }

    // Recurse into the tree to propagate info to leaves.
    for (KeySet::const_iterator key=dependencies_.begin();
         key!=dependencies_.end(); ++key) {
      S.RequireEvaluator(*key, my_tag_)->EnsureCompatibility(S);
    }
  }
}


// ---------------------------------------------------------------------------
// Updates the derivative for CompositeVectors
// ---------------------------------------------------------------------------
template<>
void EvaluatorSecondary<CompositeVector,CompositeVectorSpace>::UpdateDerivative_(State& S, const Key& wrt_key) {
  Key dmy_key = std::string("d")+my_key_+std::string("_d")+wrt_key;
  if (!S.HasData(dmy_key, my_tag_)) {
    // or create the field.  Note we have to do extra work that is normally
    // done by State in setup.
    auto& fac = S.Require<CompositeVector,CompositeVectorSpace>(dmy_key, my_tag_, my_key_);
    fac = S.Get<CompositeVector>(my_key_, my_tag_).Map();
    auto dmy_ptr = fac.Create();
    dmy_ptr->PutScalarMasterAndGhosted(0.);
    S.SetPtr(dmy_key, my_tag_, my_key_, dmy_ptr);
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_initialized();
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_io_vis(false);
    S.GetRecordW(dmy_key,my_tag_,my_key_).set_io_checkpoint(false);
  }

  auto& dmy = S.GetW<CompositeVector>(dmy_key, my_tag_, my_key_);
  dmy.PutScalarMasterAndGhosted(0.);

  // dF/dx = sum_(deps) partial F/ partial dep * ddep/dx + partial F/partial x
  for (KeySet::const_iterator dep=dependencies_.begin();
       dep!=dependencies_.end(); ++dep) {
    
    if (wrt_key == *dep) {
      // partial F / partial x
      CompositeVector tmp(dmy);
      EvaluatePartialDerivative_(S, wrt_key, tmp);
      dmy.Update(1., tmp, 1.);

    } else if (S.GetEvaluator(*dep, my_tag_)->IsDependency(S, wrt_key)) {
      // partial F / partial dep * ddep/dx
      // -- ddep/dx
      Key ddep_key = std::string("d")+*dep+std::string("_d")+wrt_key;
      auto& ddep = S.Get<CompositeVector>(ddep_key);

      // -- partial F / partial dep
      CompositeVector tmp(dmy);
      EvaluatePartialDerivative_(S, *dep, tmp);
      dmy.Multiply(1., ddep, tmp, 1.);
    }
  }
}
