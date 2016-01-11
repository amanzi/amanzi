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

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "macropore_surface_flux_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
MacroporeSurfaceFluxModel::MacroporeSurfaceFluxModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
MacroporeSurfaceFluxModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  gamma_ = plist.get<double>("gamma [-]");
  delta_ = plist.get<double>("delta [m]");
  patm_ = plist.get<double>("atmospheric pressure [Pa]", 101325.);
}


// main method
double
MacroporeSurfaceFluxModel::MacroporeSurfaceFlux(double pM, double ps, double krs, double krM, double K) const
{
  double C = gamma_ / delta_ * K * ( (pM > ps) ? krM : 0.);
  return C * (ps - pM);
}

double
MacroporeSurfaceFluxModel::DMacroporeSurfaceFluxDMacroporePressure(double pM, double ps, double krs, double krM, double K) const
{
  double C = gamma_ / delta_ * K * ( (pM > ps) ? krM : 0.);
  return -C;
}

double
MacroporeSurfaceFluxModel::DMacroporeSurfaceFluxDSurfacePressure(double pM, double ps, double krs, double krM, double K) const
{
  double C = gamma_ / delta_ * K * ( (pM > ps) ? krM : 0.);
  return C;
}

double
MacroporeSurfaceFluxModel::DMacroporeSurfaceFluxDSurfaceRelativePermeability(double pM, double ps, double krs, double krM, double K) const
{
  return 0.;
}

double
MacroporeSurfaceFluxModel::DMacroporeSurfaceFluxDMacroporeRelativePermeability(double pM, double ps, double krs, double krM, double K) const
{
  if (pM > ps)
    return gamma_ / delta_ * K * (ps - pM);
  return 0.;
}

double
MacroporeSurfaceFluxModel::DMacroporeSurfaceFluxDMacroporeAbsolutePermeability(double pM, double ps, double krs, double krM, double K) const
{
  double C = gamma_ / delta_ * ( (pM > ps) ? krM : 0.);
  return C * (ps - pM);
}

} //namespace
} //namespace
} //namespace
  
