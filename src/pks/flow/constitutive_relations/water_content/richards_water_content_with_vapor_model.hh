/*
  The richards water content with vapor model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Richards water content with vapor.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RICHARDS_WATER_CONTENT_WITH_VAPOR_MODEL_HH_
#define AMANZI_FLOW_RICHARDS_WATER_CONTENT_WITH_VAPOR_MODEL_HH_

namespace Amanzi {
namespace Flow {
namespace Relations {

class RichardsWaterContentWithVaporModel {

 public:
  explicit
  RichardsWaterContentWithVaporModel(Teuchos::ParameterList& plist);

  double WaterContent(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;

  double DWaterContentDPorosity(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  double DWaterContentDSaturationLiquid(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  double DWaterContentDSaturationGas(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  double DWaterContentDMolarDensityLiquid(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  double DWaterContentDMolarDensityGas(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  double DWaterContentDMolFractionVaporInGas(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  double DWaterContentDCellVolume(double phi, double sl, double sg, double nl, double ng, double omega, double cv) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:



};

} //namespace
} //namespace
} //namespace

#endif