/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity, given a small compressibility of rock.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_MODEL_HH_
#define AMANZI_FLOWRELATIONS_COMPRESSIBLE_POROSITY_MODEL_HH_

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class CompressiblePorosityModel {
 public:
  explicit
  CompressiblePorosityModel(Teuchos::ParameterList& plist) :
      plist_(plist) {
    InitializeFromPlist_();
  }

  virtual double Porosity(double base_poro, double pres, double patm) {
    double poro = base_poro;
    double p_over = pres - patm;
    if (p_over > cutoff_) {
      poro = base_poro + compressibility_ * ( cutoff_ / 2. + (p_over - cutoff_));
    } else if (p_over > 0.) {
      poro = base_poro + compressibility_ * (std::pow(p_over,2.) / 2. / cutoff_);
    }

    return poro > 1. ? 1. : poro;
  }

  virtual double DPorosityDPressure(double base_poro, double pres, double patm) {
    double p_over = pres - patm;
    double poro = Porosity(base_poro, pres, patm);
    if (poro == 1.) {
      return 0.;
    } else if (p_over > cutoff_) {
      return compressibility_;
    } else if (p_over > 0.) {
      return compressibility_ * p_over / cutoff_;
    }

    return 0.;
  }

  virtual double DPorosityDBasePorosity(double base_poro, double pres, double patm) {
    return pres > patm ? (Porosity(base_poro, pres, patm) > 1.0 ? 0. : 1.) : 1.;
  }

 protected:
  void InitializeFromPlist_() {
    compressibility_ = plist_.get<double>("pore compressibility");
    cutoff_ = plist_.get<double>("pore compressibility inflection point", 1000.);
  }

 protected:

  Teuchos::ParameterList plist_;
  double compressibility_;
  double cutoff_;

};

} // namespace
} // namespace
} // namespace

#endif
