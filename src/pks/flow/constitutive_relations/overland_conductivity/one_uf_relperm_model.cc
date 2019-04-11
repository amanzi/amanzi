/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the Kr associated with the unfrozen fraction of water.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "boost/math/constants/constants.hpp"
#include <cmath>

#include "dbc.hh"
#include "errors.hh"
#include "one_uf_relperm_model.hh"

namespace Amanzi {
namespace Flow {

OneUFRelPermModel::OneUFRelPermModel(Teuchos::ParameterList& plist) :
    plist_(plist),
    pi_(boost::math::constants::pi<double>()) {

  alpha_ = plist_.get<int>("unfrozen rel perm alpha", 4);
  if (alpha_ % 2 != 0) {
    Errors::Message message("Unfrozen Fraction Rel Perm: alpha must be an even integer");
    Exceptions::amanzi_throw(message);
  }

  h_cutoff_up_ = plist_.get<double>("unfrozen rel perm cutoff pressure", 10.);
  h_cutoff_dn_ = plist_.get<double>("unfrozen rel perm cutoff pressure, below", 0.);
}

double
OneUFRelPermModel::SurfaceRelPerm(double uf, double h) {
  double kr;

  if (h >= 101325. + h_cutoff_up_) {
    kr = std::pow(std::sin(pi_ * uf / 2.), alpha_);
  } else if (h <= 101325 + h_cutoff_dn_) {
    kr = 1.;
  } else {
    double fac = (h - (101325 + h_cutoff_dn_)) / (h_cutoff_up_ - h_cutoff_dn_);
    kr = (1-fac) + fac*std::pow(std::sin(pi_ * uf / 2.), alpha_);
  }  
  return kr;
}



} // namespace
} // namespace

