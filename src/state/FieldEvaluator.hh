/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/
//! A FieldEvaluator is a node in the dependency graph.
/*!

Evaluators combine data (variables) with functions (physics) providing the
values inside of that data.  Evaluators are the core concept of a model in
Amanzi-ATS: much if not all physics lives inside of an evaluator.

All evaluators are one of three types: primary, secondary, or independent
variable evaluators.  Primary variables are those that have no dependencies,
and are solved for.  Independent variables have no dependencies as well, but
are data provided by the user.  They may be constant or variable in time and
space, but they are assumed given.  Secondary variables are functions of time
and space, and of other variables.

From a software perspective, evaluators are nodes in a dependency graph, and
therefore provide the task ordering needed to evaluate the model.  Primary
variable evaluators do next to no work, but provide a target for solver to tell
the dependency graph when primary variables have been updated through the
solution/iteration of the nonlinear solve.  Independent and secondary variable
evaluators, on the other hand, implement the required processes needed to
update their data.

All field evaluators have keys they depend upon (unless they are a leaf) and
keys they provide.

All evaluator lists must provide an evaluator type, which is one of the types
registered with the evaluator factory.

.. _field-evaluator-spec:
.. admonition:: field-evaluator-spec

   * `"field evaluator type`" ``[string]`` Type registered in evaluator factory.
   * `"write checkpoint`" ``[bool]`` **true** Write this data when checkpointing.
   * `"write vis`" ``[bool]`` **true** Write this data when visualizing.

*/

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
