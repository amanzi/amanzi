/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Multiphase PK

  Secondary variable field evaluator computes product of fields
  or inverse of fields:

    eval = (f1 * f2 * ... * fn) / (g1 * g2 * ... * gm)
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
  if (my_keys_.size() == 0) {
    Key key = plist.get<std::string>("my key");
    Tag tag = make_tag(plist.get<std::string>("tag"));
    my_keys_.push_back(make_pair(key, tag));
  }
  Key domain = Keys::getDomain(my_keys_[0].first);

  if (plist_.isParameter("evaluator dependencies")) {
    Errors::Message msg;
    msg << "EvaluatorMultiplicativeReciprocal: \"" << my_keys_[0].first
        << "\" must have separate (optional) lists for multiplicative and reciprocal dependencies.";
    Exceptions::amanzi_throw(msg);
  }

  if (plist_.isParameter("multiplicative dependencies")) {
    // since dependensies is a map, we need separate maps for numerator and denominator
    const auto& names = plist_.get<Teuchos::Array<std::string>>("multiplicative dependencies");
    for (const auto& name : names) {
      Key full_name = Keys::getKey(domain, name);
      dependencies_.insert(std::make_pair(full_name, Tags::DEFAULT));
      list0_.push_back(full_name);
    }
  }

  if (plist_.isParameter("reciprocal dependencies")) {
    const auto& names = plist_.get<Teuchos::Array<std::string>>("reciprocal dependencies");
    for (const auto& name : names) {
      Key full_name = Keys::getKey(domain, name);
      dependencies_.insert(std::make_pair(full_name, Tags::DEFAULT));
      list1_.push_back(full_name);
    }
  }

  if (list0_.size() + list1_.size() == 0) {
    Errors::Message msg;
    msg << "EvaluatorMultiplicativeReciprocal for: \"" << my_keys_[0].first
        << "\" has no dependencies.";
    Exceptions::amanzi_throw(msg);
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  enforce_positivity_ = plist_.get<bool>("enforce positivity", false);
}


EvaluatorMultiplicativeReciprocal::EvaluatorMultiplicativeReciprocal(
  const EvaluatorMultiplicativeReciprocal& other)
  : EvaluatorSecondaryMonotype(other){};


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

} // namespace Amanzi
