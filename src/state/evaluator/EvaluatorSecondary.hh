/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly, see COPYRIGHT.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

//! EvaluatorSecondary is a generic function of other evaluators, potentially at
//! other tags.

/*!


*/

#ifndef STATE_EVALUATOR_SECONDARY_HH_
#define STATE_EVALUATOR_SECONDARY_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"

#include "Evaluator.hh"
#include "State.hh"

namespace Amanzi {

template <typename Data_t, typename DataFactory_t = NullFactory>
class EvaluatorSecondary : public Evaluator {

public:
  explicit EvaluatorSecondary(Teuchos::ParameterList &plist);

  EvaluatorSecondary(const EvaluatorSecondary &other) = default;

  EvaluatorSecondary &operator=(const EvaluatorSecondary &other);
  virtual Evaluator &operator=(const Evaluator &other) override;

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key reqest.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool Update(State &S, const Key &request) override final;

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State &S, const Key &request,
                                const Key &wrt_key,
                                const Key &wrt_tag) override final;

  virtual bool IsDependency(const State &S, const Key &key,
                            const Key &wrt_tag) const override final;
  virtual bool ProvidesKey(const Key &key,
                           const Key &wrt_tag) const override final;
  virtual bool IsDifferentiableWRT(const State &S, const Key &wrt_key,
                                   const Key &wrt_tag) const override final {
    return IsDependency(S, wrt_key, wrt_tag);
  }

  virtual void EnsureCompatibility(State &S) override;

  virtual std::string WriteToString() const override;

protected:
  // These do the actual work
  virtual void Evaluate_(const State &S, Data_t &result) = 0;
  virtual void EvaluatePartialDerivative_(const State &S, const Key &wrt_key,
                                          const Key &wrt_tag,
                                          Data_t &result) = 0;

  // calls Evaluate with the correct data
  virtual void Update_(State &S);

  // For general data types, we do not know how to differentiate.
  // Specializations can provide this.
  virtual void UpdateDerivative_(State &S, const Key &wrt_key,
                                 const Key &wrt_tag);
  virtual void CheckDerivative_(State &S, const Key &wrt_key,
                                const Key &wrt_tag);

protected:
  Key my_key_;
  Key my_tag_;

  KeySet requests_;
  KeyTripleSet deriv_requests_;
  KeyPairVector dependencies_;
  bool check_derivative_;

  VerboseObject vo_;
  Teuchos::ParameterList plist_;

}; // class EvaluatorSecondary

#include "EvaluatorSecondary_Impl.hh"

} // namespace Amanzi

#endif
