/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_

#include "Teuchos_ParameterList.hpp"
#include "overland_conductivity_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class ManningConductivityModel : public OverlandConductivityModel {
public:
  explicit
  ManningConductivityModel(Teuchos::ParameterList& plist);

  virtual double Conductivity(double depth, double slope, double coef);

protected:
  Teuchos::ParameterList plist_;

  double slope_regularization_;
  double manning_exp_;
  double manning_coef_;

};

} // namespace
} // namespace
} // namespace

#endif
