/*
  The ideal gas equation of state model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
    modelInitializeParamsList =   cv_ = plist.get<double>("heat capacity");
    namespaceCaps = GENERAL
    evalNameString = ideal gas equation of state
    myMethodArgs = temp_v[0][i], pres_v[0][i]
    myMethodDeclarationArgs = double temp, double pres
    myKey = density
    myKeyFirst = density
    namespace = General
    evalClassName = EosIdealGas
    paramDeclarationList =   double cv_;
    modelMethodDeclaration =   double Density(double temp, double pres) const;
    evalNameCaps = EOS_IDEAL_GAS
    myKeyMethod = Density
    evalName = eos_ideal_gas
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_ideal_gas_model.hh"

namespace Amanzi {
namespace General {
namespace Relations {

// Constructor from ParameterList
EosIdealGasModel::EosIdealGasModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
EosIdealGasModel::InitializeFromPlist_()
{
  cv_ = plist.get<double>("heat capacity");
}


// main method
double
EosIdealGasModel::Density(double temp, double pres)
{
  ASSERT(0);
}

double
EosIdealGasModel::DDensityDTemperature(double temp, double pres)
{
  ASSERT(0);
}

double
EosIdealGasModel::DDensityDPressure(double temp, double pres)
{
  ASSERT(0);
}

} //namespace
} //namespace
} //namespace
  