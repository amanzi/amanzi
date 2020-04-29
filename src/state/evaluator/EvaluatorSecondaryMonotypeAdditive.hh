/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! EvaluatorSecondaryMonotypeAdditive adds its dependencies.

/*!


*/

#ifndef STATE_EVALUATOR_ALGEBRAIC_ADDITIVE_HH_
#define STATE_EVALUATOR_ALGEBRAIC_ADDITIVE_HH_

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
class EvaluatorSecondaryMonotypeAdditive
  : public EvaluatorSecondaryMonotype<Data_t, DataFactory_t> {

 public:
  EvaluatorSecondaryMonotypeAdditive(Teuchos::ParameterList& plist);
  EvaluatorSecondaryMonotypeAdditive(const EvaluatorSecondaryMonotypeAdditive& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual std::string name() const override { return "additive"; }
  
 protected:
  virtual void
  Evaluate_(const State& S, const std::vector<Data_t*>& results) override;

  virtual void
  EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                             const Key& wrt_tag,
                             const std::vector<Data_t*>& results) override;

 protected:
  std::map<Key, double> coefs_;
  using EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::dependencies_;
  using EvaluatorSecondaryMonotype<Data_t,DataFactory_t>::my_keys_;

 private:
  static Utils::RegisteredFactory<
    Evaluator, EvaluatorSecondaryMonotypeAdditive<Data_t, DataFactory_t>>
    fac_;
};


template <typename Data_t, typename DataFactory_t>
EvaluatorSecondaryMonotypeAdditive<Data_t,DataFactory_t>::EvaluatorSecondaryMonotypeAdditive(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<Data_t, DataFactory_t>(plist)
{
  AMANZI_ASSERT(my_keys_.size() == 1);

  if (dependencies_.size() == 0) {
    Errors::Message msg("EvaluatorSecondaryMonotypeMultiplicative: empty or nonexistent \"dependencies\" list.");
    throw(msg);
  }

  for (const auto& dep : dependencies_) {
    Key coef_name = dep.first + ":" + dep.second;
    Key pname = dep.first + " coefficient";
    Key pname_full = coef_name + " coefficient";
    if (plist.isParameter(pname_full))
      coefs_[coef_name] = plist.get<double>(pname_full);
    else if (plist.isParameter(pname))
      coefs_[coef_name] = plist.get<double>(pname);
    else
      coefs_[coef_name] = 1.0;
  }
}


template <typename Data_t, typename DataFactory_t>
Teuchos::RCP<Evaluator>
EvaluatorSecondaryMonotypeAdditive<Data_t,DataFactory_t>::Clone() const
{
  return Teuchos::rcp(new EvaluatorSecondaryMonotypeAdditive(*this));
}

template <typename Data_t, typename DataFactory_t>
void
EvaluatorSecondaryMonotypeAdditive<Data_t,DataFactory_t>::Evaluate_(const State& S, const std::vector<Data_t*>& results)
{
  results[0]->putScalar(0.);
  for (const auto& dep : dependencies_) {
    const auto& term = S.Get<Data_t>(dep.first, dep.second);
    std::string coef_name = dep.first + ":" + dep.second;
    double coef = coefs_.at(coef_name);
    results[0]->update(coef, term, 1.0);
  }
}

template <typename Data_t, typename DataFactory_t>
void
EvaluatorSecondaryMonotypeAdditive<Data_t,DataFactory_t>::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key,
        const Key& wrt_tag,
        const std::vector<Data_t*>& results)
{
  Key pname_full = wrt_key + ":" + wrt_tag;
  results[0]->putScalar(coefs_[pname_full]);
}



} // namespace Amanzi

#endif
