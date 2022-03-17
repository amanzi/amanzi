/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-202x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
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
class EvaluatorSecondaryMonotypeMultiplicative : public EvaluatorSecondaryMonotype<Data_t, DataFactory_t> {
public:
  EvaluatorSecondaryMonotypeMultiplicative(Teuchos::ParameterList& plist) :
      EvaluatorSecondaryMonotype<Data_t,DataFactory_t>(plist)
  {
    coef_ = this->plist_.template get<double>("coefficient", 1.0);

    // if true, the last dependency is "divided by"
    reciprocal_ = this->plist_.template get<bool>("reciprocal", false);
  }
    

  EvaluatorSecondaryMonotypeMultiplicative(const EvaluatorSecondaryMonotypeMultiplicative& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorSecondaryMonotypeMultiplicative(*this));
  }

protected:

  virtual void Evaluate_(const State& S, const std::vector<Data_t*>& results) override {
    AMANZI_ASSERT(results.size() == 1);
    int i=0;
    results[0]->PutScalar(coef_);
    for (const auto& dep : this->dependencies_) {
      const auto& term = S.Get<Data_t>(dep.first, dep.second);
      if (reciprocal_ && this->dependencies_.size() - 1 == i) {
        results[0]->ReciprocalMultiply(1., term, *results[0], 0.);
      } else {
        results[0]->Multiply(1., term, *results[0], 0.);
      }
      ++i;
    }
  }
  
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Key& wrt_tag, const std::vector<Data_t*>& results) override {
    int i=0;
    results[0]->PutScalar(coef_);
    for (const auto& dep : this->dependencies_) {
      if (dep.first != wrt_key || dep.second != wrt_tag) {
        // not WRT
        const auto& term = S.Get<Data_t>(dep.first, dep.second);
        if (reciprocal_ && this->dependencies_.size() - 1 == i) {
          results[0]->ReciprocalMultiply(1., term, *results[0], 0.);
        } else {
          results[0]->Multiply(1., term, *results[0], 0.);
        }
      } else {
        // IS WRT
        if (reciprocal_ && this->dependencies_.size() - 1 == i) {
          //  - term ^ -2
          const auto& term = S.Get<Data_t>(dep.first, dep.second);
          results[0]->ReciprocalMultiply(-1., term, *results[0], 0.);
          results[0]->ReciprocalMultiply(1., term, *results[0], 0.);
        }
      }
      ++i;
    }
  }

 protected:
  double coef_;
  bool reciprocal_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeMultiplicative<Data_t,DataFactory_t>> fac_;
};

} // namespace Amanzi

#endif
