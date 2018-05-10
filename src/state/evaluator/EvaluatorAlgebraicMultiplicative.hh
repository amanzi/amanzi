/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! EvaluatorAlgebraicMultiplicative adds its dependencies.

/*!

  
*/

#ifndef STATE_EVALUATOR_ALGEBRAIC_MULTIPLICATIVE_HH_
#define STATE_EVALUATOR_ALGEBRAIC_MULTIPLICATIVE_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"

#include "EvaluatorAlgebraic.hh"

namespace Amanzi {

// By default, this class adds nothing on top of EvaluatorSecondary.
// Specializations can do useful things though.
template <typename Data_t, typename DataFactory_t = NullFactory>
class EvaluatorAlgebraicMultiplicative : public EvaluatorAlgebraic<Data_t, DataFactory_t> {
public:
  EvaluatorAlgebraicMultiplicative(Teuchos::ParameterList &plist) :
      EvaluatorAlgebraic<Data_t,DataFactory_t>(plist)
  {
    coef_ = this->plist_.template get<double>("coefficient", 1.0);
  }
    

  EvaluatorAlgebraicMultiplicative(const EvaluatorAlgebraicMultiplicative &other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorAlgebraicMultiplicative(*this));
  }

protected:

  virtual void Evaluate_(const State &S, Data_t &result) override {
    result.PutScalar(coef_);
    for (const auto& dep : this->dependencies_) {
      const auto& term = S.Get<Data_t>(dep.first, dep.second);
      result.Multiply(1., result, term, 0.);
    }
  }
  
  virtual void EvaluatePartialDerivative_(const State &S,
          const Key &wrt_key, const Key &wrt_tag, Data_t &result) override {
    result.PutScalar(coef_);
    for (const auto& dep : this->dependencies_) {
      if (dep.first != wrt_key || dep.second != wrt_tag) {
        const auto& term = S.Get<Data_t>(dep.first, dep.second);
        result.Multiply(1., result, term, 0.);
      }
    }
  }

 protected:
  double coef_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorAlgebraicMultiplicative<Data_t,DataFactory_t>> fac_;
  
  
};


} // namespace Amanzi

#endif
