/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  The internal energy model evaluator simply calls the IEM with 
  the correct arguments.
*/

#include "IEM_WaterVaporEvaluator.hh"

namespace Amanzi {
namespace Energy {

IEM_WaterVaporEvaluator::IEM_WaterVaporEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  // defaults work fine, this sublist need not exist
  Teuchos::ParameterList sublist = plist.sublist("IEM parameters");
  iem_ = Teuchos::rcp(new IEM_WaterVapor(sublist));

  InitializeFromPlist_();
}


IEM_WaterVaporEvaluator::IEM_WaterVaporEvaluator(
    Teuchos::ParameterList& plist, const Teuchos::RCP<IEM_WaterVapor>& iem)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      iem_(iem) {
  InitializeFromPlist_();
}


IEM_WaterVaporEvaluator::IEM_WaterVaporEvaluator(const IEM_WaterVaporEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      temp_key_(other.temp_key_),
      mol_frac_key_(other.mol_frac_key_),
      iem_(other.iem_) {};


Teuchos::RCP<Evaluator>
IEM_WaterVaporEvaluator::Clone() const {
  return Teuchos::rcp(new IEM_WaterVaporEvaluator(*this));
}


void IEM_WaterVaporEvaluator::InitializeFromPlist_()
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(make_pair(plist_.get<std::string>("internal energy key"), Tags::DEFAULT));
  }

  // Set up my dependencies.
  std::string domain = Keys::getDomain(my_keys_[0].first);
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(std::make_pair(temp_key_, Tags::DEFAULT));

  // -- molar fraction of water vapor in the gaseous phase
  mol_frac_key_ = plist_.get<std::string>("vapor molar fraction key", Keys::getKey(domain, "molar_fraction_gas"));
  dependencies_.insert(std::make_pair(mol_frac_key_, Tags::DEFAULT));
}


void IEM_WaterVaporEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  auto temp = S.GetPtr<CompositeVector>(temp_key_);
  auto mol_frac = S.GetPtr<CompositeVector>(mol_frac_key_);

  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
    const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

    int ncomp = results[0]->size(*comp);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i], molfrac_v[0][i]);
    }
  }
}


void IEM_WaterVaporEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  auto temp = S.GetPtr<CompositeVector>(temp_key_);
  auto mol_frac = S.GetPtr<CompositeVector>(mol_frac_key_);

  if (wrt_key == temp_key_) {
    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

      int ncomp = results[0]->size(*comp);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else if (wrt_key == mol_frac_key_) {
    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

      int ncomp = results[0]->size(*comp);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDomega(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace Energy
}  // namespace Amanzi
