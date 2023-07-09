/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#ifndef AMANZI_TABULAR_MODEL_HH_
#define AMANZI_TABULAR_MODEL_HH_

#include <memory>

#include <boost/math/tools/roots.hpp>
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "SplinedCurve.hh"

#include "WRM.hh"

namespace Amanzi {
namespace Flow {

class WRM_tabular : public WRM {
  struct Tol {
    Tol(double eps) : eps_(eps){};
    bool operator()(double a, double b) const { return std::abs(a - b) <= eps_; }
    double eps_;
  };

  struct F {
    F(double s, const WRM_tabular* wrm) : s_(s), wrm_(wrm){};
    double operator()(double pc) const { return wrm_->saturation(pc) - s_; }
    double s_;
    const WRM_tabular* wrm_;
  };

 public:
  explicit WRM_tabular(Teuchos::ParameterList& plist);
  ~WRM_tabular(){};

  // required methods from the base class
  double k_relative(double pc) const { return spline_kr_->Value(pc); }
  double saturation(double pc) const { return spline_sat_->Value(pc); }
  double dSdPc(double pc) const { return spline_sat_->Derivative(pc); }
  double capillaryPressure(double saturation) const;
  double residualSaturation() const { return sr_; }
  double dKdPc(double pc) const { return spline_kr_->Derivative(pc); }

  // access 
  boost::uintmax_t get_itrs() { return itrs_; }

 private:
  std::shared_ptr<Utils::SplinedCurve> spline_kr_, spline_sat_;
  double pc0_, sr_;
  mutable boost::uintmax_t itrs_;

  static Utils::RegisteredFactory<WRM, WRM_tabular> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
