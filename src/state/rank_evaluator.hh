/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for MPI ranks.

------------------------------------------------------------------------- */

#ifndef AMANZI_COMM_RANK_FIELD_EVALUATOR_
#define AMANZI_COMM_RANK_FIELD_EVALUATOR_

#include "field_evaluator.hh"
#include "field_evaluator_factory.hh"

namespace Amanzi {

class RankModel : public FieldEvaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  RankModel(Teuchos::ParameterList& plist);

  RankModel(const RankModel& other);

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
  std::string my_key_;
  std::string my_mesh_;
  Teuchos::ParameterList plist_;
  bool computed_once_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,RankModel> fac_;
};

} // namespace


#endif
