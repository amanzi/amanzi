/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation, slope, and updates pres + elev.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "snow_skin_potential_evaluator.hh"

namespace Amanzi {
namespace Flow {

SnowSkinPotentialEvaluator::SnowSkinPotentialEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  if (my_key_.empty()) 
    my_key_ = plist_.get<std::string>("potential key", "snow_skin_potential");

  pd_key_ = plist_.get<std::string>("ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);
  sd_key_ = plist_.get<std::string>("snow depth key", "snow_depth");
  dependencies_.insert(sd_key_);
  precip_key_ = plist_.get<std::string>("precipitation snow key", "precipitation_snow");
  dependencies_.insert(precip_key_);
  elev_key_ = plist_.get<std::string>("elevation key", "elevation");
  dependencies_.insert(elev_key_);

  factor_ = plist_.get<double>("dt factor");
}


SnowSkinPotentialEvaluator::SnowSkinPotentialEvaluator(const SnowSkinPotentialEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    elev_key_(other.elev_key_),
    pd_key_(other.pd_key_),
    sd_key_(other.sd_key_),
    precip_key_(other.precip_key_),
    factor_(other.factor_)
{};


Teuchos::RCP<FieldEvaluator>
SnowSkinPotentialEvaluator::Clone() const {
  return Teuchos::rcp(new SnowSkinPotentialEvaluator(*this));
}


void SnowSkinPotentialEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  // update pressure + elevation
  Teuchos::RCP<const CompositeVector> pd = S->GetFieldData(pd_key_);
  Teuchos::RCP<const CompositeVector> sd = S->GetFieldData(sd_key_);
  Teuchos::RCP<const CompositeVector> precip = S->GetFieldData(precip_key_);
  Teuchos::RCP<const CompositeVector> elev = S->GetFieldData(elev_key_);

  
  result->Update(1.0, *elev, 1.0, *pd, 0.0);
  result->Update(1.0, *sd, factor_, *precip, 1.0);
}


// This is hopefully never called?
void SnowSkinPotentialEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  ASSERT(0);
  result->PutScalar(1.0);
}

} //namespace
} //namespace
