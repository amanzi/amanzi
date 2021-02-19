/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A field evaluator with no dependencies specified by a function.

------------------------------------------------------------------------- */

#ifndef AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFILE_
#define AMANZI_INDEPENDENT_FIELD_EVALUATOR_FROMFILE_

#include "independent_variable_field_evaluator.hh"
#include "FieldEvaluator_Factory.hh"

namespace Amanzi {

class IndependentVariableFieldEvaluatorFromFile :
      public IndependentVariableFieldEvaluator
{

 public:

  // ---------------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------------
  explicit
  IndependentVariableFieldEvaluatorFromFile(Teuchos::ParameterList& plist);
  IndependentVariableFieldEvaluatorFromFile(const IndependentVariableFieldEvaluatorFromFile& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  // ---------------------------------------------------------------------------
  // Update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void UpdateField_(const Teuchos::Ptr<State>& S);

  void LoadFile_(int i);

  void Interpolate_(double time, const Teuchos::Ptr<CompositeVector>& v);

 protected:
  double t_before_, t_after_;
  Teuchos::RCP<CompositeVector> val_before_, val_after_;
  std::string filename_;
  std::vector<double> times_;
  int current_interval_;

  std::string meshname_;
  std::string compname_;
  std::string varname_;
  std::string locname_;
  int ndofs_;

  Teuchos::RCP<Function> time_func_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,IndependentVariableFieldEvaluatorFromFile> fac_;
};

} // namespace


#endif
