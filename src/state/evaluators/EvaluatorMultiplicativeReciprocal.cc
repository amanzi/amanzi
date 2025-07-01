/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <utility>

#include "Evaluator.hh"
#include "EvaluatorMultiplicativeReciprocal.hh"

namespace Amanzi {

/* ******************************************************************
* Two constructors.
****************************************************************** */
EvaluatorMultiplicativeReciprocal::EvaluatorMultiplicativeReciprocal(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  if (plist_.isParameter("evaluator dependencies")) {
    Errors::Message msg;
    msg << "EvaluatorMultiplicativeReciprocal: \"" << my_keys_.front().first
        << "\" must have separate (optional) lists for multiplicative and reciprocal dependencies.";
    Exceptions::amanzi_throw(msg);
  }
  if (plist_.isParameter("multiplicative dependencies") ||
      plist_.isParameter("reciprocal dependencies")) {
    Errors::Message msg;
    msg << "EvaluatorMultiplicativeReciprocal: \"" << my_keys_.front().first
        << "\" no longer accepts option \"multiplicative dependencies\" or \"reciprocal dependencies\""
        << "-- please use \"multiplicative dependency keys\" or \"multiplicative dependency key suffixes\""
        << " (respectively reciprocal) instead.";
  }

  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  const Teuchos::Array<std::string> empty_array;
  {
    list0_ = Keys::readKeys(plist_, domain, "multiplicative dependency", &empty_array);
    for (const auto& key : list0_) dependencies_.insert({key, tag});
  }
  {
    list1_ = Keys::readKeys(plist_, domain, "reciprocal dependency", &empty_array);
    for (const auto& key : list1_) dependencies_.insert({key, tag});
  }

  if (list0_.size() + list1_.size() == 0) {
    Errors::Message msg;
    msg << "EvaluatorMultiplicativeReciprocal for: \"" << my_keys_.front().first
        << "\" has no dependencies.";
    Exceptions::amanzi_throw(msg);
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  enforce_positivity_ = plist_.get<bool>("enforce positivity", false);
}


EvaluatorMultiplicativeReciprocal::EvaluatorMultiplicativeReciprocal(
  const EvaluatorMultiplicativeReciprocal& other)
  : EvaluatorSecondaryMonotype(other)
{};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator>
EvaluatorMultiplicativeReciprocal::Clone() const
{
  return Teuchos::rcp(new EvaluatorMultiplicativeReciprocal(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
EvaluatorMultiplicativeReciprocal::Evaluate_(const State& S,
                                             const std::vector<CompositeVector*>& results)
{
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    auto& result_c = *results[0]->ViewComponent(*comp);
    int ndofs = result_c.MyLength();

    result_c.PutScalar(coef_);

    for (auto it = list0_.begin(); it != list0_.end(); ++it) {
      const auto& factor_c = *S.Get<CompositeVector>(*it, Tags::DEFAULT).ViewComponent(*comp);
      for (int c = 0; c != ndofs; ++c) result_c[0][c] *= factor_c[0][c];
    }

    for (auto it = list1_.begin(); it != list1_.end(); ++it) {
      const auto& factor_c = *S.Get<CompositeVector>(*it, Tags::DEFAULT).ViewComponent(*comp);
      for (int c = 0; c != ndofs; ++c) result_c[0][c] /= factor_c[0][c];
    }

    if (enforce_positivity_) {
      for (int c = 0; c != ndofs; ++c) { result_c[0][c] = std::max(result_c[0][c], 0.0); }
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void
EvaluatorMultiplicativeReciprocal::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    auto& result_c = *results[0]->ViewComponent(*comp);
    int ncells = result_c.MyLength();

    result_c.PutScalar(1.0);
    for (auto it = list0_.begin(); it != list0_.end(); ++it) {
      const auto& factor_c = *S.Get<CompositeVector>(*it, Tags::DEFAULT).ViewComponent(*comp);
      if (*it != wrt_key)
        for (int c = 0; c != ncells; ++c) result_c[0][c] *= factor_c[0][c];
    }

    for (auto it = list1_.begin(); it != list1_.end(); ++it) {
      const auto& factor_c = *S.Get<CompositeVector>(*it, Tags::DEFAULT).ViewComponent(*comp);
      if (*it == wrt_key)
        for (int c = 0; c != ncells; ++c) result_c[0][c] /= -factor_c[0][c] * factor_c[0][c];
      else
        for (int c = 0; c != ncells; ++c) result_c[0][c] /= factor_c[0][c];
    }
  }
}


/* ******************************************************************
* Units are calculated if field has none. Otherwise, units are compared.
****************************************************************** */
void
EvaluatorMultiplicativeReciprocal::EnsureCompatibility_Units_(State& S)
{
  std::string data("-");
  Utils::Units system;
  auto& r = S.GetRecordSetW(my_keys_[0].first);

  for (auto it = list0_.begin(); it != list0_.end(); ++it) {
    auto tmp = S.GetRecordSet(*it).units();
    data = system.MultiplyUnits(data, tmp);
  }
  for (auto it = list1_.begin(); it != list1_.end(); ++it) {
    auto tmp = S.GetRecordSet(*it).units();
    data = system.DivideUnits(data, tmp);
  }

  auto tmp = r.units();
  if (tmp != "")
    AMANZI_ASSERT(system.CompareUnits(tmp, data));
  else
    r.set_units(data);
}

} // namespace Amanzi
