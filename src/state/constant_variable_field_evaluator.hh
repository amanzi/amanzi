/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  A field evaluator with no dependencies initialized from file and have no 
  time dependencies.
*/

#ifndef AMANZI_INDEPENDENT_CONSTANT_FIELD_EVALUATOR_
#define AMANZI_INDEPENDENT_CONSTANT_FIELD_EVALUATOR_

#include "CompositeVectorFunction.hh"
#include "FieldEvaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class ConstantVariableFieldEvaluator : public FieldEvaluator {

public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  ConstantVariableFieldEvaluator(Teuchos::ParameterList& plist);
  ConstantVariableFieldEvaluator(const ConstantVariableFieldEvaluator& other);

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

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual std::string WriteToString() const;

 protected:
  Key my_key_;

  KeySet requests_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,ConstantVariableFieldEvaluator> fac_;
};

}  // namespace Amanzi

#endif
