/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Arcos

License: BSD
Author: Ethan Coon

An evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_EVALUATOR_PRIMARY_
#define AMANZI_STATE_EVALUATOR_PRIMARY_

#include <memory>

#include "Evaluator.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorPrimary : public Evaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  EvaluatorPrimary(Teuchos::ParameterList& plist);
  EvaluatorPrimary(const EvaluatorPrimary& other) = default;

  // virtual destructor
  virtual ~EvaluatorPrimary() {}

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  EvaluatorPrimary& operator=(const EvaluatorPrimary& other);
  virtual Evaluator& operator=(const Evaluator& other) override;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of the evaluator.
  //
  // Updates the data, if needed.  Returns true if the value of the data has
  // changed since the last request for an update.
  // ---------------------------------------------------------------------------
  virtual bool Update(State& S, const Key& request) override;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of derivatives of evaluator.
  //
  // Updates the derivative, if needed.  Returns true if the value of the
  // derivative with respect to wrt_key has changed since the last request for
  // an update.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State& S,
          const Key& request, const Key& wrt_key) override;

  virtual bool IsDependency(const State& S, const Key& key) const override;
  virtual bool ProvidesKey(const Key& key) const override;

  // ---------------------------------------------------------------------------
  // How a PK informs this leaf of the tree that it has changed.
  //
  // Effectively this simply tosses the request history, so that the next
  // requests will say this has changed.
  // ---------------------------------------------------------------------------
  virtual void SetChanged();

  virtual void EnsureCompatibility(State& S) override {}
  virtual std::string WriteToString() const override; 

 protected:
  Key my_key_;
  Key my_tag_;
  KeySet requests_;
  KeyPairSet deriv_requests_;

  VerboseObject vo_;
  
 private:
  static Utils::RegisteredFactory<Evaluator,EvaluatorPrimary> fac_;

};

} // namespace

#endif
