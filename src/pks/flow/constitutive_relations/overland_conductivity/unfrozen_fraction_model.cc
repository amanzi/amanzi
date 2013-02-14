/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen fraction of water.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "unfrozen_fraction_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

UnfrozenFractionModel::UnfrozenFractionModel(Teuchos::ParameterList& plist) :
    plist_(plist) {
  halfwidth_ = plist_.get<double>("transition width", 0.1) / 2.;
  T0_ = plist_.get<double>("freezing point", 273.15);
  pi_ = boost::math::constants::pi<double>();
}

double UnfrozenFractionModel::UnfrozenFraction(double temp) {
  double adj_temp = temp - T0_;
  double uf;
  if (adj_temp > 0.) {
    uf = 1.;
  } else if (adj_temp < -2*halfwidth_) {
    uf = 0.;
  } else {
    uf = (std::sin(pi_/2. * (adj_temp + halfwidth_)/halfwidth_) + 1.)/2.;
  }
  return uf;
}

double UnfrozenFractionModel::DUnfrozenFractionDT(double temp) {
  double adj_temp = temp - T0_;
  double duf;
  if (adj_temp > 0.) {
    duf = 0.;
  } else if (adj_temp < -2*halfwidth_) {
    duf = 0.;
  } else {
    duf = std::cos(pi_/2. * (adj_temp + halfwidth_)/halfwidth_)/2. * pi_/2. / halfwidth_;
  }
  return duf;
}


} // namespace
} // namespace
} // namespace
