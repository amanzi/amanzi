/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Wraps a PDE_Accumulation to be an Evaluator.

/*!

Lots of options here, document me!  
  
*/

#include "Evaluator_PDE_Accumulation.hh"

namespace Amanzi {

Evaluator_PDE_Accumulation::Evaluator_PDE_Accumulation(Teuchos::ParameterList &plist)
    : EvaluatorAlgebraic<CompositeVector,CompositeVectorSpace>(plist)
{
  conserved_key_ = plist.get<std::string>("conserved quantity");
  tag_old_ = plist.get<std::string>("tag old");
  tag_new_ = plist.get<std::string>("tag new");

  dependencies_.emplace_back(std::make_pair(conserved_key_, tag_old_));
  dependencies_.emplace_back(std::make_pair(conserved_key_, tag_new_));
}

void
Evaluator_PDE_Accumulation::EnsureCompatibility(State &S) {
  EvaluatorAlgebraic<CompositeVector,CompositeVectorSpace>::EnsureCompatibility(S);
  S.Require<double>("time", tag_old_);
  S.Require<double>("time", tag_new_);
}
  
void
Evaluator_PDE_Accumulation::Evaluate_(const State &S, CompositeVector &result) {
  double dt = S.Get<double>("time", tag_new_) - S.Get<double>("time", tag_old_);
  result.Update(1./dt, S.Get<CompositeVector>(conserved_key_, tag_new_),
                -1./dt, S.Get<CompositeVector>(conserved_key_, tag_old_), 0.);
}

void
Evaluator_PDE_Accumulation::EvaluatePartialDerivative_(const State &S,
        const Key &wrt_key, const Key &wrt_tag, CompositeVector &result) {
  ASSERT(wrt_key == conserved_key_);
  double dt = S.Get<double>("time", tag_new_) - S.Get<double>("time", tag_old_);
  if (wrt_tag == tag_old_) {
    result.PutScalar(-1./dt);
  } else if (wrt_tag == tag_new_) {
    result.PutScalar(1./dt);
  } else {
    ASSERT(0);
  }
}

}  // namespace Amanzi

