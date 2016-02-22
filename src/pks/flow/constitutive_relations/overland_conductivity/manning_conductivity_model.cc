/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

ManningConductivityModel::ManningConductivityModel(Teuchos::ParameterList& plist) :
    plist_(plist) {

  slope_regularization_ = plist_.get<double>("slope regularization epsilon", 1.e-8);
  manning_exp_ = plist_.get<double>("Manning exponent");
}

double ManningConductivityModel::Conductivity(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;

  double exponent = manning_exp_ + 1.0;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return std::pow(std::max(depth,0.), exponent) / scaling;
}

double ManningConductivityModel::DConductivityDDepth(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;
  double exponent = manning_exp_ + 1.0;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return std::pow(std::max(depth,0.), exponent - 1.) * exponent / scaling;
}


} // namespace
} // namespace
} // namespace
