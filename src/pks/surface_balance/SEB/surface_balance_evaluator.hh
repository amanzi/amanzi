/*
  SurfaceBalanceEvaluator evaluates SEB as a nonlinear source instead of as a PK.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_EVALUATOR_HH_

#include "factory.hh"
#include "secondary_variables_field_evaluator.hh"

namespace Amanzi {

class Debugger;

namespace SurfaceBalance {

class SurfaceBalanceEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  SurfaceBalanceEvaluator(Teuchos::ParameterList& plist);

  SurfaceBalanceEvaluator(const SurfaceBalanceEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results);

 protected:
  double min_wind_speed_;
  double snow_ground_trans_;
  double albedo_trans_;
  bool initialized_;

  Teuchos::RCP<Debugger> db_;


 private:
  static Utils::RegisteredFactory<FieldEvaluator,SurfaceBalanceEvaluator> reg_;
};

} // namespace
} // namespace

#endif
