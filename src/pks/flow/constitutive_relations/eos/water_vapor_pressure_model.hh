/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_WATER_VAPOR_PRESSURE_HH_
#define _FLOWRELATIONS_EOS_WATER_VAPOR_PRESSURE_HH_


#include "Teuchos_ParameterList.hpp"
#include "factory.hh"
#include "vapor_pressure_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WaterVaporPressureModel : public VaporPressureModel {

public:
  explicit WaterVaporPressureModel(Teuchos::ParameterList& plist);

  virtual double SaturatedVaporPressure(double T);
  virtual double DSaturatedVaporPressureDT(double T);

private:
  Teuchos::ParameterList plist_;
  const double ka0_;
  const double ka_, kb_, kc_, kd_;

  static Utils::RegisteredFactory<VaporPressureModel,WaterVaporPressureModel> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
