/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies solved for by a PK.

------------------------------------------------------------------------- */

#ifndef AMANZI_PRIMARY_VARIABLE_FIELD_EVALUATOR_
#define AMANZI_PRIMARY_VARIABLE_FIELD_EVALUATOR_

#include "FieldEvaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class PrimaryVariableFieldEvaluator : public FieldEvaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  PrimaryVariableFieldEvaluator(Teuchos::ParameterList& plist);
  PrimaryVariableFieldEvaluator(const PrimaryVariableFieldEvaluator& other);

  // virtual destructor
  virtual ~PrimaryVariableFieldEvaluator() {}

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;
  virtual void operator=(const FieldEvaluator& other);

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key reqest.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
                                         Key request, Key wrt_key);

  virtual bool IsDependency(const Teuchos::Ptr<State>& S, Key key) const;
  virtual bool ProvidesKey(Key key) const;


  // ---------------------------------------------------------------------------
  // How a PK informs this leaf of the tree that it has changed.
  //
  // Effectively this simply tosses the request history, so that the next
  // requests will say this has changed.
  // ---------------------------------------------------------------------------
  virtual void SetFieldAsChanged(const Teuchos::Ptr<State>& S);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) {}

  virtual std::string WriteToString() const;

  // special flags
  virtual bool IsPrimary() { return true; }

 protected:
  Key my_key_;
  KeySet requests_;
  KeyPairSet deriv_requests_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,PrimaryVariableFieldEvaluator> fac_;

};

} // namespace

#endif
