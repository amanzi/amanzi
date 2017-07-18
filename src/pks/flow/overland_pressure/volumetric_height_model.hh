/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates Volumetric height(pressure) for subgrid model

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOW_VOLUMETRIC_HEIGHT_MODEL_
#define AMANZI_FLOW_VOLUMETRIC_HEIGHT_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class VolumetricHeightModel {
public:
  explicit
  VolumetricHeightModel(Teuchos::ParameterList& plist) : plist_(plist) {}

  double Height(double pd, double delta_max, double delta_ex) {
    double d =  std::pow(pd,2)*(2*delta_max - 3*delta_ex)/std::pow(delta_max,2) + std::pow(pd,3)*(2*delta_ex - delta_max)/std::pow(delta_max,3);
    return d;
  }
  
  double DHeightDPressure(double pd, double delta_max, double delta_ex) {
    std::cout<<"DH_DP for the subgrid model is not implemented yet!! \n";abort();
    return 1. / (1.);
  }

  double DHeightDRho(double pres, double rho, double p_atm, double g_z) {
    std::cout<<"DH_DRho for the subgrid model is not implemented yet!! \n";abort();
    return 0.;//-(pres - p_atm) / (rho * rho * g_z);
  }

protected:
  Teuchos::ParameterList plist_;

};

} // namespace
} // namespace

#endif
