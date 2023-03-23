/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Energy

  Interface for a thermal conductivity model with one phase.
*/

#include "dbc.hh"
#include "EOSFactory.hh"
#include "EOS_ThermalConductivity.hh"
#include "EOS_Utils.hh"

#include "TCMEvaluator_OnePhase.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
TCMEvaluator_OnePhase::TCMEvaluator_OnePhase(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(
      make_pair(plist_.get<std::string>("thermal conductivity key"), Tags::DEFAULT));
  }
  auto prefix = Keys::getDomainPrefix(my_keys_[0].first);

  temperature_key_ = plist_.get<std::string>("temperature key", prefix + "temperature");
  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));

  porosity_key_ = plist_.get<std::string>("porosity key", prefix + "porosity");
  dependencies_.insert(std::make_pair(porosity_key_, Tags::DEFAULT));

  AmanziEOS::EOSFactory<AmanziEOS::EOS_ThermalConductivity> eos_fac;
  auto& sublist = plist_.sublist("thermal conductivity parameters");
  tc_ = eos_fac.Create(sublist);

  k_rock_ = sublist.get<double>("thermal conductivity of rock");
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
TCMEvaluator_OnePhase::TCMEvaluator_OnePhase(const TCMEvaluator_OnePhase& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    tc_(other.tc_),
    temperature_key_(other.temperature_key_){};


/* ******************************************************************
* TBW.
****************************************************************** */
Teuchos::RCP<Evaluator>
TCMEvaluator_OnePhase::Clone() const
{
  return Teuchos::rcp(new TCMEvaluator_OnePhase(*this));
}


/* ******************************************************************
* Evaluator body.
****************************************************************** */
void
TCMEvaluator_OnePhase::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // pull out the dependencies
  const auto& temp_c = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");
  const auto& poro_c = *S.Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  Epetra_MultiVector& result_c = *results[0]->ViewComponent("cell");

  int ierr(0);
  int ncomp = results[0]->size("cell", false);
  for (int i = 0; i != ncomp; ++i) {
    double k_liq = tc_->ThermalConductivity(temp_c[0][i], 0.0); // no models with pressure yet
    ierr = std::max(ierr, tc_->error_code());

    double phi = poro_c[0][i];
    result_c[0][i] = phi * k_liq + (1.0 - phi) * k_rock_;
  }
  AmanziEOS::ErrorAnalysis(S.Get<CompositeVector>(temperature_key_).Comm(), ierr, tc_->error_msg());
}


/* ******************************************************************
* Evaluator of derivarives.
****************************************************************** */
void
TCMEvaluator_OnePhase::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& results)
{
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const auto& temp_v = *S.Get<CompositeVector>(temperature_key_).ViewComponent(*comp);
    const auto& poro_v = *S.Get<CompositeVector>(porosity_key_).ViewComponent(*comp);
    Epetra_MultiVector& result_v = *results[0]->ViewComponent(*comp);
    int ncells = results[0]->size(*comp);

    if (wrt_key == porosity_key_) {
      for (int i = 0; i != ncells; ++i) {
        double k_liq = tc_->ThermalConductivity(temp_v[0][i], 0.0);
        result_v[0][i] = k_liq - k_rock_;
      }
    } else if (wrt_key == temperature_key_) {
      for (int i = 0; i != ncells; ++i) {
        double k_liq = tc_->DThermalConductivityDT(temp_v[0][i], 0.0);
        result_v[0][i] = poro_v[0][i] * k_liq;
      }
    }
  }
}

} // namespace Energy
} // namespace Amanzi
