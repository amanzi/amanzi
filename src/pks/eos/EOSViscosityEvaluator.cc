/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  ViscosityEvaluator is the interface between state/data and the 
  model, a VPM.
*/

#include "EOS_Utils.hh"
#include "EOSFactory.hh"
#include "EOSViscosityEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

EOSViscosityEvaluator::EOSViscosityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  // my keys
  if (my_keys_.size() == 0)
    my_keys_.push_back(
      std::make_pair(plist_.get<std::string>("viscosity key", "viscosity_liquid"), Tags::DEFAULT));

  // set up my dependencies
  Key domain = Keys::getDomain(my_keys_[0].first);
  tag_ = Tags::DEFAULT;

  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(std::make_pair(temp_key_, tag_));
  pres_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(std::make_pair(pres_key_, tag_));

  // Construct my Viscosity model
  AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory<EOS_Viscosity> eos_fac;
  visc_ = eos_fac.Create(plist_.sublist("EOS parameters"));
}


EOSViscosityEvaluator::EOSViscosityEvaluator(const EOSViscosityEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    visc_(other.visc_),
    temp_key_(other.temp_key_){};


Teuchos::RCP<Evaluator>
EOSViscosityEvaluator::Clone() const
{
  return Teuchos::rcp(new EOSViscosityEvaluator(*this));
}


void
EOSViscosityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  auto temp = S.GetPtr<CompositeVector>(temp_key_, tag_);
  auto pres = S.GetPtr<CompositeVector>(pres_key_, tag_);

  // evaluate p_s / p_atm
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int ierr = 0;
    int count = results[0]->size(*comp);
    for (int i = 0; i != count; ++i) {
      result_v[0][i] = visc_->Viscosity(temp_v[0][i], pres_v[0][i]);
      ierr = std::max(ierr, visc_->error_code());
    }
    ErrorAnalysis(temp->Comm(), ierr, visc_->error_msg());
  }
}


void
EOSViscosityEvaluator::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& results)
{
  auto temp = S.GetPtr<CompositeVector>(temp_key_, tag_);
  auto pres = S.GetPtr<CompositeVector>(pres_key_, tag_);

  // evaluate d/dT( p_s / p_atm )
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp));
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp));
    Epetra_MultiVector& result_v = *(results[0]->ViewComponent(*comp));

    int count = results[0]->size(*comp);
    for (int i = 0; i != count; ++i) {
      result_v[0][i] = visc_->DViscosityDT(temp_v[0][i], pres_v[0][i]);
    }
  }
}

} // namespace AmanziEOS
} // namespace Amanzi
