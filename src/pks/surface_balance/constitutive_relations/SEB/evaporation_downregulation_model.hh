/*
  The evaporation downregulation via soil resistance model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Downregulates evaporation from a potential.

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SEB_EVAPORATION_DOWNREGULATION_MODEL_HH_
#define AMANZI_SEB_EVAPORATION_DOWNREGULATION_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class EvaporationDownregulationModel {

 public:
  explicit
  EvaporationDownregulationModel(Teuchos::ParameterList& plist);

  double Evaporation(double sg, double poro, double pot_evap) const;

  double DEvaporationDSaturationGas(double sg, double poro, double pot_evap) const;
  double DEvaporationDPorosity(double sg, double poro, double pot_evap) const;
  double DEvaporationDPotentialEvaporation(double sg, double poro, double pot_evap) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double dess_dz_;
  double Clapp_Horn_b_;

};

} //namespace
} //namespace
} //namespace

#endif