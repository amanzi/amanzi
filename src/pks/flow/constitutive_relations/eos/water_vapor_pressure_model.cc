/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include <cmath>
#include "water_vapor_pressure_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<VaporPressureModel,WaterVaporPressureModel> WaterVaporPressureModel::factory_("water vapor over water/ice");

WaterVaporPressureModel::WaterVaporPressureModel(Teuchos::ParameterList& plist) :
  plist_(plist),
  ka0_(16.635764),
  ka_(-6096.9385),
  kb_(-2.7111933e-2),
  kc_(1.673952e-5),
  kd_(2.433502) {}

double WaterVaporPressureModel::SaturatedVaporPressure(double T) {
  return 100.0*exp(ka0_ + ka_/T + (kb_ + kc_*T)*T + kd_*log(T));
};

double WaterVaporPressureModel::DSaturatedVaporPressureDT(double T) {
  return 100.0*exp(ka0_ + ka_/T + (kb_ + kc_*T)*T + kd_*log(T))
    * (-ka_/(T*T) + kb_ + kc_*T + kd_/T);
};


} //namespace
} //namespace
} //namespace
