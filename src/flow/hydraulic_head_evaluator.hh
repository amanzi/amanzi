#ifndef AMANZI_HYDRAULIC_HEAD_EVALUATOR_
#define AMANZI_HYDRAULIC_HEAD_EVALUATOR_

#include "factory.hh"
#include "FieldEvaluator.hh"
#include "FieldEvaluator_Factory.hh"


namespace Amanzi {
namespace Flow {

class HydraulicHeadEvaluator : public FieldEvaluator {

public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  HydraulicHeadEvaluator(Teuchos::ParameterList& plist);
  HydraulicHeadEvaluator(const HydraulicHeadEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;
  virtual void operator=(const FieldEvaluator& other) {};

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

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S);

 protected:
  Key my_key_;
  Key my_mesh_;
  bool computed_once_;
  bool communicate_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,HydraulicHeadEvaluator> fac_;
};


} //namespace
} //namespace


#endif
