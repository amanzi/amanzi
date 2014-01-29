/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the Kr associated with the unfrozen fraction of water.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "boost/math/constants/constants.hpp"
#include <cmath>

#include "dbc.hh"
#include "errors.hh"
#include "unfrozen_fraction_relperm_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

UnfrozenFractionRelPermModel::UnfrozenFractionRelPermModel(Teuchos::ParameterList& plist) :
    plist_(plist),
    pi_(boost::math::constants::pi<double>()) {
  alpha_ = plist_.get<int>("unfrozen rel perm alpha", 4);
  if (alpha_ % 2 != 0) {
    Errors::Message message("Unfrozen Fraction Rel Perm: alpha must be an even integer");
    Exceptions::amanzi_throw(message);
  }

  min_val_ = plist_.get<double>("minimum rel perm cutoff", 1.e-30);
  h_cutoff_ = plist_.get<double>("unfrozen rel perm cutoff height", 0.01);
}

double
UnfrozenFractionRelPermModel::UnfrozenFractionRelPerm(double uf, double h) {
  double kr = std::pow(std::sin(pi_ * uf / 2.), alpha_);
  if (h < h_cutoff_) {
    // double fac = h / h_cutoff_;
    // kr = kr * fac + (1-fac); // kr --> 1  as h --> 0

    double fac = std::pow(std::sin(pi_ * (h/h_cutoff_) / 2.), 2);
    kr *= fac; // kr --> 0 as h --> 0

  }
  return std::max(kr, min_val_);
}

double
UnfrozenFractionRelPermModel::DUnfrozenFractionRelPermDUnfrozenFraction(double uf, double h) {
  ASSERT(0);
  return 0.;
}


} // namespace
} // namespace
} // namespace
