/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <cmath>
#include <string>

#include "Brent.hh"
#include "errors.hh"
#include "FlowDefs.hh"
#include "WRM_tabular.hh"

namespace Amanzi {
namespace Flow {


/* ******************************************************************
* Constructors
****************************************************************** */
WRM_tabular::WRM_tabular(Teuchos::ParameterList& plist) : tol_(1.0e-11)
{
  auto pc = plist.get<Teuchos::Array<double>>("cap pressure").toVector();
  auto kr = plist.get<Teuchos::Array<double>>("permeability").toVector();
  auto sat = plist.get<Teuchos::Array<double>>("saturation").toVector();

  tol_ = plist.get<double>("tolerance", 1.0e-11);

  if (pc.size() != kr.size() || pc.size() != sat.size()) {
    Errors::Message msg;
    msg << "tabular: tabulated arrays have different size.";
    Exceptions::amanzi_throw(msg);
  }

  spline_kr_ = std::make_shared<Utils::SplinedCurve>(
    pc,
    kr,
    std::make_pair(Utils::SplinedCurve::SplineExtrapolation_t::CONSTANT,
                   Utils::SplinedCurve::SplineExtrapolation_t::CONSTANT),
    true);

  spline_sat_ = std::make_shared<Utils::SplinedCurve>(
    pc,
    sat,
    std::make_pair(Utils::SplinedCurve::SplineExtrapolation_t::CONSTANT,
                   Utils::SplinedCurve::SplineExtrapolation_t::CONSTANT),
    true);

  pc0_ = pc[pc.size() - 1];
  sr_ = sat[sat.size() - 1];
}


/* ******************************************************************
* Capillary pressure as a function of saturation.
****************************************************************** */
double
WRM_tabular::capillaryPressure(double s) const
{
  F f(s, this);

  itrs_ = 100;
  return Utils::computeRootBrent(f, 0.0, pc0_, tol_, &itrs_);
}

} // namespace Flow
} // namespace Amanzi
