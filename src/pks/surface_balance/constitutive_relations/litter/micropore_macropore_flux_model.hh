/*
  The micropore-macropore flux model is an algebraic model with dependencies.

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

#ifndef AMANZI_SURFACEBALANCE_MICROPORE_MACROPORE_FLUX_MODEL_HH_
#define AMANZI_SURFACEBALANCE_MICROPORE_MACROPORE_FLUX_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MicroporeMacroporeFluxModel {

 public:
  explicit
  MicroporeMacroporeFluxModel(Teuchos::ParameterList& plist);

  double MicroporeMacroporeFlux(double pm, double pM, double krM, double krm, double K) const;

  double DMicroporeMacroporeFluxDMicroporePressure(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDPressure(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDRelativePermeability(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDMicroporeRelativePermeability(double pm, double pM, double krM, double krm, double K) const;
  double DMicroporeMacroporeFluxDMicroporeAbsolutePermeability(double pm, double pM, double krM, double krm, double K) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double gamma_;
  double delta_;

};

} //namespace
} //namespace
} //namespace

#endif
