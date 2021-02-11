/*
  The plant wilting factor model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Wilting factor.

Beta, or the water availability factor, or the plant wilting factor.

Beta =  (p_closed - p) / (p_closed - p_open)

where p is the capillary pressure or soil mafic potential, and closed
and open indicate the values at which stomates are fully open or fully
closed (the wilting point).


    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_PLANT_WILTING_FACTOR_MODEL_HH_
#define AMANZI_FLOW_PLANT_WILTING_FACTOR_MODEL_HH_

namespace Amanzi {
namespace LandCover {
namespace Relations {

class PlantWiltingFactorModel {

 public:
  explicit
  PlantWiltingFactorModel(Teuchos::ParameterList& plist);

  double PlantWiltingFactor(double pc) const;

  double DPlantWiltingFactorDCapillaryPressureGasLiq(double pc) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double pc_o_;
  double pc_c_;

};

} //namespace
} //namespace
} //namespace

#endif