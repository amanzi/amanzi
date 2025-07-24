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
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), n_dofs_(-1)
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
        << "\" no longer accepts option \"multiplicative dependencies\" or \"reciprocal "
           "dependencies\""
        << "-- please use \"multiplicative dependency keys\" or \"multiplicative dependency key "
           "suffixes\""
        << " (respectively reciprocal) instead.";
  }

  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  const Teuchos::Array<std::string> empty_array;
  {
    list0_ = Keys::readKeys(plist_, domain, "multiplicative dependency", &empty_array);
    for (const auto& key : list0_) dependencies_.insert({ key, tag });
  }
  {
    list1_ = Keys::readKeys(plist_, domain, "reciprocal dependency", &empty_array);
    for (const auto& key : list1_) dependencies_.insert({ key, tag });
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
  : EvaluatorSecondaryMonotype(other) {};


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
  for (const auto& comp : *results[0]) {
    int n_dofs = results[0]->NumVectors(comp);

    auto& result_c = *results[0]->ViewComponent(comp);
    result_c.PutScalar(coef_);

    for (const auto& it : list0_) {
      const auto& factor_c = *S.Get<CompositeVector>(it, Tags::DEFAULT).ViewComponent(comp);
      if (factor_c.NumVectors() == n_dofs) {
        for (int i = 0; i != n_dofs; ++i) {
          for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] *= factor_c[i][c];
        }
      } else {
        for (int i = 0; i != n_dofs; ++i) {
          for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] *= factor_c[0][c];
        }
      }
    }

    for (const auto& it : list1_) {
      const auto& factor_c = *S.Get<CompositeVector>(it, Tags::DEFAULT).ViewComponent(comp);
      if (factor_c.NumVectors() == n_dofs) {
        for (int i = 0; i != n_dofs; ++i) {
          for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] /= factor_c[i][c];
        }
      } else {
        for (int i = 0; i != n_dofs; ++i) {
          for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] /= factor_c[0][c];
        }
      }
    }

    if (enforce_positivity_) {
      for (int c = 0; c != result_c.MyLength(); ++c) {
        result_c[0][c] = std::max(result_c[0][c], 0.0);
      }
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
  for (const auto& comp : *results[0]) {
    auto& result_c = *results[0]->ViewComponent(comp);
    result_c.PutScalar(coef_);

    for (const auto& it : list0_) {
      const auto& factor_c = *S.Get<CompositeVector>(it, Tags::DEFAULT).ViewComponent(comp);
      if (it != wrt_key) {
        if (factor_c.NumVectors() == n_dofs_) {
          for (int i = 0; i != n_dofs_; ++i) {
            for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] *= factor_c[i][c];
          }
        } else {
          for (int i = 0; i != n_dofs_; ++i) {
            for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] *= factor_c[0][c];
          }
        }
      }
    }

    for (const auto& it : list1_) {
      const auto& factor_c = *S.Get<CompositeVector>(it, Tags::DEFAULT).ViewComponent(comp);
      if (it == wrt_key) {
        if (factor_c.NumVectors() == n_dofs_) {
          for (int i = 0; i != n_dofs_; ++i) {
            for (int c = 0; c != result_c.MyLength(); ++c)
              result_c[i][c] /= -factor_c[i][c] * factor_c[i][c];
          }
        } else {
          for (int i = 0; i != n_dofs_; ++i) {
            for (int c = 0; c != result_c.MyLength(); ++c)
              result_c[i][c] /= -factor_c[0][c] * factor_c[0][c];
          }
        }

      } else {
        if (factor_c.NumVectors() == n_dofs_) {
          for (int i = 0; i != n_dofs_; ++i) {
            for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] /= factor_c[i][c];
          }
        } else {
          for (int i = 0; i != n_dofs_; ++i) {
            for (int c = 0; c != result_c.MyLength() ; ++c) result_c[i][c] /= factor_c[0][c];
          }
        }
      }
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

  bool valid_units = true;
  for (auto it = list0_.begin(); it != list0_.end(); ++it) {
    auto tmp = S.GetRecordSet(*it).units();
    if (tmp != "") {
      data = system.MultiplyUnits(data, tmp);
    } else {
      valid_units = false;
      break;
    }
  }
  for (auto it = list1_.begin(); it != list1_.end(); ++it) {
    auto tmp = S.GetRecordSet(*it).units();
    if (tmp != "") {
      data = system.DivideUnits(data, tmp);
    } else {
      valid_units = false;
      break;
    }
  }

  if (valid_units) {
    auto& r = S.GetRecordSetW(my_keys_[0].first);
    auto tmp = r.units();
    if (tmp != "") AMANZI_ASSERT(system.CompareUnits(tmp, data));
    else r.set_units(data);
  }
}

void
EvaluatorMultiplicativeReciprocal::EnsureCompatibility_ToDeps_(State& S)
{
  auto akeytag = my_keys_.front();
  const auto& my_fac =
    S.Require<CompositeVector, CompositeVectorSpace>(akeytag.first, akeytag.second);

  if (my_fac.Mesh() != Teuchos::null) {
    for (const auto& key_tag : dependencies_) {
      auto& dep_fac =
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second);
      dep_fac.SetMesh(my_fac.Mesh());

      for (const auto& comp : my_fac) {
        AmanziMesh::Entity_kind loc = my_fac.Location(comp);
        int n_dofs = my_fac.NumVectors(comp);

        if (n_dofs_ > 0) {
          AMANZI_ASSERT(n_dofs_ == n_dofs);
        } else {
          n_dofs_ = n_dofs;
        }

        if (dep_fac.HasComponent(comp)) {
          int dep_n_dofs = dep_fac.NumVectors(comp);
          if (dep_n_dofs == 1 || dep_n_dofs == n_dofs_) {
            dep_fac.AddComponent(comp, loc, dep_n_dofs);
          } else {
            AMANZI_ASSERT(false);
          }
        }
      }
    }
  }
}


} // namespace Amanzi
