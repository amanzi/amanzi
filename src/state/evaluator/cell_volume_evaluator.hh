/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator for an unchanging cell volume.

------------------------------------------------------------------------- */

#ifndef AMANZI_CELL_VOLUME_FIELD_EVALUATOR_
#define AMANZI_CELL_VOLUME_FIELD_EVALUATOR_

#include "Evaluator.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class CellVolumeEvaluator : public Evaluator {

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  CellVolumeEvaluator(Teuchos::ParameterList& plist);
  CellVolumeEvaluator(const CellVolumeEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;
  virtual void operator=(const Evaluator& other) {};

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
  virtual bool UpdateDerivative(const Teuchos::Ptr<State>& S,
                                         Key request, Key wrt_key);

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
  Key my_key_;
  Key my_mesh_;
  bool computed_once_;
  bool communicate_;

 private:
  static Utils::RegisteredFactory<Evaluator,CellVolumeEvaluator> fac_;
};

} // namespace


#endif
