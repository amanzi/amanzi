/*
  The richards water content with vapor model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content with vapor.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "richards_water_content_with_vapor_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
RichardsWaterContentWithVaporModel::RichardsWaterContentWithVaporModel(Teuchos::ParameterList& plist)
{
  InitializeFromPlist_(plist);
}


// Initialize parameters
void
RichardsWaterContentWithVaporModel::InitializeFromPlist_(Teuchos::ParameterList& plist)
{

}


// main method
double
RichardsWaterContentWithVaporModel::WaterContent(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*phi*(ng*omega*sg + nl*sl);
}

double
RichardsWaterContentWithVaporModel::DWaterContentDPorosity(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*(ng*omega*sg + nl*sl);
}

double
RichardsWaterContentWithVaporModel::DWaterContentDSaturationLiquid(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*nl*phi;
}

double
RichardsWaterContentWithVaporModel::DWaterContentDSaturationGas(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*ng*omega*phi;
}

double
RichardsWaterContentWithVaporModel::DWaterContentDMolarDensityLiquid(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*phi*sl;
}

double
RichardsWaterContentWithVaporModel::DWaterContentDMolarDensityGas(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*omega*phi*sg;
}

double
RichardsWaterContentWithVaporModel::DWaterContentDMolFractionVaporInGas(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return cv*ng*phi*sg;
}

double
RichardsWaterContentWithVaporModel::DWaterContentDCellVolume(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const
{
  return phi*(ng*omega*sg + nl*sl);
}

} //namespace
} //namespace
} //namespace
  