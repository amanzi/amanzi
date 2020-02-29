/*
  The micropore-macropore flux evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    evalName = micropore_macropore_flux
    modelMethodDeclaration =   double MicroporeMacroporeFlux(double pm, double pM, double krM, double krm, double K) const;
    namespaceCaps = SURFACEBALANCE
    namespace = SurfaceBalance
    evalNameCaps = MICROPORE_MACROPORE_FLUX
    myMethodArgs = pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]
    myKeyMethod = MicroporeMacroporeFlux
    myKeyFirst = micropore
    evalNameString = micropore-macropore flux
    myMethodDeclarationArgs = double pm, double pM, double krM, double krm, double K
    evalClassName = MicroporeMacroporeFlux
    myKey = micropore_macropore_flux
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_MICROPORE_MACROPORE_FLUX_EVALUATOR_HH_
#define AMANZI_SURFACEBALANCE_MICROPORE_MACROPORE_FLUX_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MicroporeMacroporeFluxModel;

class MicroporeMacroporeFluxEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  MicroporeMacroporeFluxEvaluator(Teuchos::ParameterList& plist);
  MicroporeMacroporeFluxEvaluator(const MicroporeMacroporeFluxEvaluator& other);

  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  Teuchos::RCP<MicroporeMacroporeFluxModel> get_model() { return model_; }

 protected:
  void InitializeFromPlist_();

  Key pm_key_;
  Key pM_key_;
  Key krM_key_;
  Key krm_key_;
  Key K_key_;

  Teuchos::RCP<MicroporeMacroporeFluxModel> model_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,MicroporeMacroporeFluxEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
