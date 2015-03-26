/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface class for a FieldEvaluator.  A FieldEvaluator is a node in the Phalanx-like
dependency tree.

------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_FIELD_EVALUATOR_HH_
#define AMANZI_STATE_FIELD_EVALUATOR_HH_

#include <string>
#include <vector>
#include <ostream>

#include "VerboseObject.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "State.hh"

namespace Amanzi {

class FieldEvaluator {

 public:
  explicit
  FieldEvaluator(Teuchos::ParameterList& plist);

  FieldEvaluator(const FieldEvaluator& other);

  // virtual destructor
  virtual ~FieldEvaluator() {}

  virtual Teuchos::RCP<FieldEvaluator> Clone() const = 0;
  virtual void operator=(const FieldEvaluator& other) = 0;

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key request.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request) = 0;

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool HasFieldDerivativeChanged(const Teuchos::Ptr<State>& S,
          Key requester, Key wrt_key) = 0;

  // Update the field if needed.
  //  virtual UpdateField(const Teuchos::Ptr<State>& S) = 0;

  // Update the field's derivative with respect to wrt_key if needed.
  //  virtual UpdateFieldDerivative(const Teuchos::Ptr<State>& S, Key wrt_key) = 0;


  virtual bool IsDependency(const Teuchos::Ptr<State>& S, Key key) const = 0;
  virtual bool ProvidesKey(Key key) const = 0;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) = 0;

  virtual std::string WriteToString() const = 0;

  friend std::ostream& operator<<(std::ostream&, const FieldEvaluator&);

 protected:
  // parameter list for the object
  Teuchos::ParameterList plist_;

  // VerboseObject for output
  Teuchos::RCP<VerboseObject> vo_;

}; // class FieldEvaluator



} // namespace Amanzi

#endif
