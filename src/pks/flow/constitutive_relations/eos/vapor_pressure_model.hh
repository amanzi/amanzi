/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas air with a molar fraction of water vapor.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef _FLOWRELATIONS_EOS_VAPOR_PRESSURE_HH_
#define _FLOWRELATIONS_EOS_VAPOR_PRESSURE_HH_

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class VaporPressureModel {

public:
  // explicit VaporPressureModel(Teuchos::ParameterList& plist);

  virtual double SaturatedVaporPressure(double T) = 0;
  virtual double DSaturatedVaporPressureDT(double T) = 0;

};

} //namespace
} //namespace
} //namespace

#endif
