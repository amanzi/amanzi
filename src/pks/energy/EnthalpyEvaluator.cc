/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for enthalpy, H = U + p / eta.
*/

#include "EnthalpyEvaluator.hh"
#include "IEMEvaluator.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor of a secondary evaluator.
****************************************************************** */
EnthalpyEvaluator::EnthalpyEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  tag_ = Tags::DEFAULT;
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("enthalpy key"), tag_));
  }
  auto domain = Keys::getDomain(my_keys_[0].first);  // include dash

  // Set up my dependencies.
  // -- internal energy
  ie_key_ = plist_.get<std::string>("internal energy key", Keys::getKey(domain, "internal_energy_liquid"));
  dependencies_.insert(std::make_pair(ie_key_, Tags::DEFAULT));

  // -- pressure work
  include_work_ = plist_.get<bool>("include work term", true);

  if (include_work_) {
    pressure_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
    dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));

    mol_density_key_ = plist_.get<std::string>("molar density key", Keys::getKey(domain, "molar_density_liquid"));
    dependencies_.insert(std::make_pair(mol_density_key_, Tags::DEFAULT));
  }
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
EnthalpyEvaluator::EnthalpyEvaluator(const EnthalpyEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      pressure_key_(other.pressure_key_),
      mol_density_key_(other.mol_density_key_),
      ie_key_(other.ie_key_),
      include_work_(other.include_work_) {};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<Evaluator> EnthalpyEvaluator::Clone() const {
  return Teuchos::rcp(new EnthalpyEvaluator(*this));
}


/* ******************************************************************
* Evaluator's body.
****************************************************************** */
void EnthalpyEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  auto u_l = S.GetPtr<CompositeVector>(ie_key_, tag_);
  *results[0] = *u_l;

  if (include_work_) {
    auto pres = S.GetPtr<CompositeVector>(pressure_key_, tag_);
    auto n_l = S.GetPtr<CompositeVector>(mol_density_key_, tag_);

    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

      int ncomp = results[0]->size(*comp);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] += pres_v[0][i] / nl_v[0][i];
      }
    }
  }
}


/* ******************************************************************
* Evaluator for derivatives.
****************************************************************** */
void EnthalpyEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  if (wrt_key == ie_key_) {
    results[0]->PutScalar(1.0);

  } else if (wrt_key == pressure_key_) {
    AMANZI_ASSERT(include_work_);
    
    auto n_l = S.GetPtr<CompositeVector>(mol_density_key_, tag_);

    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

      int ncomp = results[0]->size(*comp);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = 1.0 / nl_v[0][i];
      }
    }

  } else if (wrt_key == mol_density_key_) {
    AMANZI_ASSERT(include_work_);
    
    auto pres = S.GetPtr<CompositeVector>(pressure_key_, tag_);
    auto n_l = S.GetPtr<CompositeVector>(mol_density_key_, tag_);

    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const Epetra_MultiVector& nl_v = *n_l->ViewComponent(*comp);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);

      int ncomp = results[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = -pres_v[0][i] / (nl_v[0][i] * nl_v[0][i]);
      }
    }
  }
}


/* ******************************************************************
* Evaluation at a point
****************************************************************** */
double EnthalpyEvaluator::EvaluateFieldSingle(
    const Teuchos::Ptr<State>& S, int c, double T, double p)
{
  double tmp = Teuchos::rcp_dynamic_cast<IEMEvaluator>(S->GetEvaluatorPtr(ie_key_, Tags::DEFAULT))->EvaluateFieldSingle(c, T, p);
  if (include_work_) {
    const auto& nl_c = *S->Get<CompositeVector>(mol_density_key_).ViewComponent("cell", true);
    tmp += p / nl_c[0][c];
  }
  return tmp;
}

}  // namespace Energy
}  // namespace Amanzi
