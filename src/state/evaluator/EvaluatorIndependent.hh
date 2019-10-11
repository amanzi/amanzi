/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/


//! Base class for evaluators with no dependencies.

/*!

Note this cannot be constructed, and instead is constructed from a derived
class such as EvaluatorIndependentFromFunction or EvaluatorIndependentFromFile.

*/

#ifndef AMANZI_INDEPENDENT_EVALUATOR_
#define AMANZI_INDEPENDENT_EVALUATOR_

#include "errors.hh"
#include "Evaluator.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

//
// Dummy class, does everything but know the type, which is required to
// EnsureCompatibility.  This is never used, instead the below templated one
// is.
//
class EvaluatorIndependent_ : public Evaluator {
 public:
  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  explicit EvaluatorIndependent_(Teuchos::ParameterList& plist);
  EvaluatorIndependent_(const EvaluatorIndependent_& other) = default;

  EvaluatorIndependent_& operator=(const EvaluatorIndependent_& other);
  virtual Evaluator& operator=(const Evaluator& other) override;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of the evaluator.
  //
  // Updates the data, if needed.  Returns true if the value of the data has
  // changed since the last request for an update.
  // ---------------------------------------------------------------------------
  virtual bool Update(State& S, const Key& request) override final;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of derivatives of evaluator.
  //
  // Updates the derivative, if needed.  Returns true if the value of the
  // derivative with respect to wrt_key has changed since the last request for
  // an update.
  // ---------------------------------------------------------------------------
  virtual bool
  UpdateDerivative(State& S, const Key& request, const Key& wrt_key,
                   const Key& wrt_tag) override final;

  virtual bool IsDependency(const State& S, const Key& key,
                            const Key& tag) const override final;
  virtual bool ProvidesKey(const Key& key, const Key& tag) const override final;
  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key,
                                   const Key& wrt_tag) const override final;

  virtual void EnsureCompatibility(State& S) override;
  //  virtual void EnsureCompatibleDerivative(State &S, const Key& wrt_key,
  //  const Key& wrt_tag) override;

  virtual std::string WriteToString() const override;

 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) = 0;

 protected:
  Key my_key_;
  Key my_tag_;

  double time_;
  bool temporally_variable_;
  bool computed_once_;

  KeySet requests_;
  Teuchos::ParameterList plist_;
  VerboseObject vo_;
};

template <class Data_t, class DataFactory_t = NullFactory>
class EvaluatorIndependent : public EvaluatorIndependent_ {
 public:
  using EvaluatorIndependent_::EvaluatorIndependent_;

  virtual void EnsureCompatibility(State& S) override
  {
    // Require the field and claim ownership.
    S.Require<Data_t, DataFactory_t>(my_key_, my_tag_, my_key_);
  }
};

} // namespace Amanzi

#endif
