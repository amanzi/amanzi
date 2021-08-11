/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! A field evaluator providing user-supplied data as a leaf node.
/*!

Independent variables are data provided by the user.  Independent variable
evaluators are therefore leaf nodes in the dependency graph (they have no
dependencies themselves) but do know how to update their data with values as a
function of space and time.

There are three ways for the user to provide that data: as static, constant
values, from an Amanzi Function_, or from file.

*/

#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_

#include "CompositeVectorFunction.hh"
#include "FieldEvaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class IndependentVariableFieldEvaluator : public FieldEvaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  IndependentVariableFieldEvaluator(Teuchos::ParameterList& plist);
  IndependentVariableFieldEvaluator(const IndependentVariableFieldEvaluator& other);

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

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  virtual std::string WriteToString() const;

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S) = 0;

 protected:
  Key my_key_;

  double time_;
  bool temporally_variable_;
  bool computed_once_;

  KeySet requests_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluator> fac_;
};

} // namespace


#endif
