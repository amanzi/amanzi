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

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "micropore_macropore_flux_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
MicroporeMacroporeFluxModel::MicroporeMacroporeFluxModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
MicroporeMacroporeFluxModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  gamma_ = plist.get<double>("gamma [-]");
  delta_ = plist.get<double>("delta [m]");
}


// main method
double
MicroporeMacroporeFluxModel::MicroporeMacroporeFlux(double pm, double pM, double krM, double krm, double K) const
{
  double C = gamma_ / delta_ * K * (pm > pM ? krm : krM);
  return C * (pM - pm);
}

double
MicroporeMacroporeFluxModel::DMicroporeMacroporeFluxDMicroporePressure(double pm, double pM, double krM, double krm, double K) const
{
  double C = gamma_ / delta_ * K * (pm > pM ? krm : krM);
  return -C;
}

double
MicroporeMacroporeFluxModel::DMicroporeMacroporeFluxDPressure(double pm, double pM, double krM, double krm, double K) const
{
  double C = gamma_ / delta_ * K * (pm > pM ? krm : krM);
  return C;
}

double
MicroporeMacroporeFluxModel::DMicroporeMacroporeFluxDRelativePermeability(double pm, double pM, double krM, double krm, double K) const
{
  double C = gamma_ / delta_ * K * (pm > pM ? 0. : 1.);
  return C * (pM - pm);
}

double
MicroporeMacroporeFluxModel::DMicroporeMacroporeFluxDMicroporeRelativePermeability(double pm, double pM, double krM, double krm, double K) const
{
  double C = gamma_ / delta_ * K * (pm > pM ? 1. : 0.);
  return C * (pM - pm);
}

double
MicroporeMacroporeFluxModel::DMicroporeMacroporeFluxDMicroporeAbsolutePermeability(double pm, double pM, double krM, double krm, double K) const
{
  double C = gamma_ / delta_ * (pm > pM ? krm : krM);
  return C * (pM - pm);
}

} //namespace
} //namespace
} //namespace
  
