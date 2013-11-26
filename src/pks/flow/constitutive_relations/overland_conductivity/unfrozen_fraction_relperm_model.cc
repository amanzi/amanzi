/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the Kr associated with the unfrozen fraction of water.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "boost/math/constants/constants.hpp"
#include <cmath>

#include "dbc.hh"
#include "unfrozen_fraction_relperm_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

UnfrozenFractionRelPermModel::UnfrozenFractionRelPermModel(Teuchos::ParameterList& plist) :
    plist_(plist),
    pi_(boost::math::constants::pi<double>()) {}

double
UnfrozenFractionRelPermModel::UnfrozenFractionRelPerm(double uf, double h) {
  return std::pow(std::sin(pi_ * uf / 2.), 4);
}

double
UnfrozenFractionRelPermModel::DUnfrozenFractionRelPermDUnfrozenFraction(double uf, double h) {
  ASSERT(0);
  return 0.;
}


} // namespace
} // namespace
} // namespace
