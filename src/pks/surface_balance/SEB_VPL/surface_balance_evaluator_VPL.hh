/*
  SurfaceBalanceEvaluator evaluates SEB as a nonlinear source instead of as a PK.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_EVALUATOR_VPL_HH_
#define AMANZI_SURFACEBALANCE_EVALUATOR_VPL_HH_

#include "factory.hh"
#include "secondary_variable_field_evaluator.hh"


namespace Amanzi {

class Debugger;

namespace SurfaceBalance {

class SurfaceBalanceEvaluatorVPL : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  SurfaceBalanceEvaluatorVPL(Teuchos::ParameterList& plist);

  SurfaceBalanceEvaluatorVPL(const SurfaceBalanceEvaluatorVPL& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

 protected:

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  double min_wind_speed_;
  double snow_ground_trans_;
  double albedo_trans_;
  bool initialized_;

  Teuchos::RCP<Debugger> db_;
  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SurfaceBalanceEvaluatorVPL> reg_;
};

} // namespace
} // namespace

#endif
