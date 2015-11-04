/*
  The macropore-surface flux model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
    evalName = macropore_surface_flux
    modelMethodDeclaration =   double MacroporeSurfaceFlux(double pM, double ps, double krs, double krM, double K) const;
    namespaceCaps = SURFACEBALANCE
    namespace = SurfaceBalance
    evalNameCaps = MACROPORE_SURFACE_FLUX
    myMethodArgs = pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]
    myKeyMethod = MacroporeSurfaceFlux
    myKeyFirst = macropore
    evalNameString = macropore-surface flux
    myMethodDeclarationArgs = double pM, double ps, double krs, double krM, double K
    evalClassName = MacroporeSurfaceFlux
    myKey = macropore_surface_flux
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_MACROPORE_SURFACE_FLUX_MODEL_HH_
#define AMANZI_SURFACEBALANCE_MACROPORE_SURFACE_FLUX_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class MacroporeSurfaceFluxModel {

 public:
  explicit
  MacroporeSurfaceFluxModel(Teuchos::ParameterList& plist);

  double MacroporeSurfaceFlux(double pM, double ps, double krs, double krM, double K) const;

  double DMacroporeSurfaceFluxDMacroporePressure(double pM, double ps, double krs, double krM, double K) const;
  double DMacroporeSurfaceFluxDSurfacePressure(double pM, double ps, double krs, double krM, double K) const;
  double DMacroporeSurfaceFluxDSurfaceRelativePermeability(double pM, double ps, double krs, double krM, double K) const;
  double DMacroporeSurfaceFluxDMacroporeRelativePermeability(double pM, double ps, double krs, double krM, double K) const;
  double DMacroporeSurfaceFluxDMacroporeAbsolutePermeability(double pM, double ps, double krs, double krM, double K) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double gamma_;
  double delta_;
  double patm_;

};

} //namespace
} //namespace
} //namespace

#endif
