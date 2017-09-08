/*
  The rooting depth fraction model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.


    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "rooting_depth_fraction_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
RootingDepthFractionModel::RootingDepthFractionModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
RootingDepthFractionModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  a_ = plist.get<double>("alpha", 7);
  b_ = plist.get<double>("beta", 1.75);
  z_max_ = plist.get<double>("max rooting depth [m]", 2.0);
}


// main method
double
RootingDepthFractionModel::RootingDepthFraction(double z) const
{
  return ((z <= z_max_) ? (
   0.5*a_*exp(-a_*z) + 0.5*b_*exp(-b_*z)
)
: (
   0.0
));
}

double
RootingDepthFractionModel::DRootingDepthFractionDDepth(double z) const
{
  return ((z <= z_max_) ? (
   -0.5*pow(a_, 2)*exp(-a_*z) - 0.5*pow(b_, 2)*exp(-b_*z)
)
: (
   0
));
}

} //namespace
} //namespace
} //namespace
  