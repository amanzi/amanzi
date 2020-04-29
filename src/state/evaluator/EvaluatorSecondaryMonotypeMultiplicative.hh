/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! EvaluatorSecondaryMonotypeMultiplicative adds its dependencies.

/*!


*/

#ifndef STATE_EVALUATOR_ALGEBRAIC_MULTIPLICATIVE_HH_
#define STATE_EVALUATOR_ALGEBRAIC_MULTIPLICATIVE_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

// By default, this class adds nothing on top of EvaluatorSecondary.
// Specializations can do useful things though.
template <typename Data_t, typename DataFactory_t = NullFactory>
class EvaluatorSecondaryMonotypeMultiplicative
  : public EvaluatorSecondaryMonotype<Data_t, DataFactory_t> {
 public:
  EvaluatorSecondaryMonotypeMultiplicative(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<Data_t, DataFactory_t>(plist)
  {
    if (dependencies_.size() == 0) {
      Errors::Message msg("EvaluatorSecondaryMonotypeMultiplicative: empty or nonexistent \"dependencies\" list.");
      throw(msg);
    }

    coef_ = plist_.template get<double>("coefficient", 1.0);

    // if true, the last dependency is "divided by"
    reciprocal_ = plist_.template get<bool>("reciprocal", false);
  }


  EvaluatorSecondaryMonotypeMultiplicative(
    const EvaluatorSecondaryMonotypeMultiplicative& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorSecondaryMonotypeMultiplicative(*this));
  }

  virtual std::string name() const override { return "multiplicative"; }
  
 protected:
  virtual void
  Evaluate_(const State& S, const std::vector<Data_t*>& results) override
  {
    AMANZI_ASSERT(results.size() == 1);
    AMANZI_ASSERT(dependencies_.size() >= 1);

    const auto& dep_last = dependencies_.back();
    const auto& term_last = S.Get<Data_t>(dep_last.first, dep_last.second);
    if (reciprocal_) {
      results[0]->reciprocal(term_last);
      if (coef_ != 1.0) results[0]->scale(coef_);
    } else {
      results[0]->update(coef_, term_last, 0.);
    }

    for (int lcv=0; lcv!=dependencies_.size()-1; ++lcv) {
      const auto& dep = dependencies_[lcv];
      const auto& term = S.Get<Data_t>(dep.first, dep.second);
      results[0]->elementWiseMultiply(1., term, *results[0], 0.);
    }
  }

  virtual void
  EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                             const Key& wrt_tag,
                             const std::vector<Data_t*>& results) override
  {
    auto wrt = std::make_pair(wrt_key, wrt_tag);
    
    const auto& dep_last = dependencies_.back();
    const auto& term_last = S.Get<Data_t>(dep_last.first, dep_last.second);
    if (reciprocal_) {
      if (dep_last != wrt) {
        results[0]->reciprocal(term_last);
        if (coef_ != 1.0) results[0]->scale(coef_);
      } else {
        results[0]->elementWiseMultiply(-1/coef_, term_last, term_last, 0.);
        results[0]->reciprocal(*results[0]);
      }
    } else {
      if (dep_last != wrt) {
        results[0]->update(coef_, term_last, 0.);
      } else {
        results[0]->putScalar(coef_);
      }
    }

    for (int lcv=0; lcv!=dependencies_.size()-1; ++lcv) {
      const auto& dep = dependencies_[lcv];
      if (dep != wrt) {
        const auto& term = S.Get<Data_t>(dep.first, dep.second);
        results[0]->elementWiseMultiply(1., term, *results[0], 0.);
      }
    }
  }

 protected:
  double coef_;
  bool reciprocal_;
  using EvaluatorSecondaryMonotype<Data_t, DataFactory_t>::plist_;
  using EvaluatorSecondaryMonotype<Data_t, DataFactory_t>::dependencies_;

 private:
  static Utils::RegisteredFactory<
    Evaluator, EvaluatorSecondaryMonotypeMultiplicative<Data_t, DataFactory_t>>
    fac_;
};


} // namespace Amanzi

#endif
