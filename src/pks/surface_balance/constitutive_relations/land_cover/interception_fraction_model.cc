/*
  The interception fraction model is an algebraic model with dependencies.

  Generated via evaluator_generator with:

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "interception_fraction_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
InterceptionFractionModel::InterceptionFractionModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
InterceptionFractionModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{
  alpha_ = plist.get<double>("leaf area interception fraction [-]", 0.25);
}


// main method
double
InterceptionFractionModel::InterceptionFraction(double ai) const
{
  return alpha_*(1 - exp(-0.5*ai));
}

double
InterceptionFractionModel::DInterceptionFractionDAreaIndex(double ai) const
{
  return 0.5*alpha_*exp(-0.5*ai);
}

} //namespace
} //namespace
} //namespace
  
