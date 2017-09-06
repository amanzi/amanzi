/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator base for secondary variables.  A Evaluator is a node in the
Phalanx-like dependency tree.

Secondary variable evaluators, such as equations of state, water retention evaluators,
internal energy evaluators, etc should inherit this class, implementing the
missing Update_() and UpdateFieldDerivative_() methods.

------------------------------------------------------------------------- */

#ifndef STATE_EVALUATOR_SECONDARY_HH_
#define STATE_EVALUATOR_SECONDARY_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "exceptions.hh"

#include "State.hh"
#include "Evaluator.hh"

namespace Amanzi {

template<typename Data_t, typename DataFactory_t=NullFactory>
class EvaluatorSecondary : public Evaluator {

 public:
  explicit
  EvaluatorSecondary(Teuchos::ParameterList& plist);

  EvaluatorSecondary(const EvaluatorSecondary& other) = default;

  EvaluatorSecondary& operator=(const EvaluatorSecondary& other);
  virtual Evaluator& operator=(const Evaluator& other) override;

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key reqest.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool Update(State& S, const Key& request) override;

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State& S,
          const Key& request, const Key& wrt_key) override;

  virtual bool IsDependency(const State& S, const Key& key) const override;
  virtual bool ProvidesKey(const Key& key) const override;

  virtual void EnsureCompatibility(State& S) override;

  virtual std::string WriteToString() const override;

 protected:
  // These do the actual work
  virtual void Evaluate_(const State& S,
                         Data_t& result) = 0;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, Data_t& result) = 0;

  // calls Evaluate with the correct data
  virtual void Update_(State& S);

  // For general data types, we do not know how to differentiate.
  // Specializations can provide this.
  virtual void UpdateDerivative_(State& S, const Key& wrt_key);

  virtual void CheckDerivative_(State& S, const Key& wrt_key);
  
 protected:
  Key my_key_;
  Key my_tag_;

  KeySet requests_;
  KeyPairSet deriv_requests_;
  KeySet dependencies_;
  bool check_derivative_;

  VerboseObject vo_;
  Teuchos::ParameterList plist_;
  
}; // class EvaluatorSecondary

#include "EvaluatorSecondary_Impl.hh"
#include "EvaluatorSecondaryDouble_Impl.hh"
#include "EvaluatorSecondaryCompositeVector_Impl.hh"

} // namespace Amanzi

#endif
