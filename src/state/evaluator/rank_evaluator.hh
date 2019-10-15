/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_COMM_RANK_FIELD_EVALUATOR_
#define AMANZI_COMM_RANK_FIELD_EVALUATOR_

#include "Evaluator.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class RankEvaluator : public Evaluator {
 public:
  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit RankEvaluator(Teuchos::ParameterList& plist);

  RankEvaluator(const RankEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;
  virtual void operator=(const Evaluator& other){};

  // ---------------------------------------------------------------------------
  // Answers the question, has this Field changed since it was last requested
  // for Field Key reqest.  Updates the field if needed.
  // ---------------------------------------------------------------------------
  virtual bool Update(const Teuchos::Ptr<State>& S, Key request);

  // ---------------------------------------------------------------------------
  // Answers the question, Has This Field's derivative with respect to Key
  // wrt_key changed since it was last requested for Field Key reqest.
  // Updates the derivative if needed.
  // ---------------------------------------------------------------------------
  virtual bool
  UpdateDerivative(const Teuchos::Ptr<State>& S, Key request, Key wrt_key);

  virtual bool IsDependency(const Teuchos::Ptr<State>& S, Key key) const;
  virtual bool ProvidesKey(Key key) const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

  virtual std::string WriteToString() const;

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(const Teuchos::Ptr<State>& S);

 protected:
  std::string my_key_;
  std::string my_mesh_;
  Teuchos::ParameterList plist_;
  bool computed_once_;

 private:
  static Utils::RegisteredFactory<Evaluator, RankEvaluator> fac_;
};

} // namespace Amanzi

#endif
