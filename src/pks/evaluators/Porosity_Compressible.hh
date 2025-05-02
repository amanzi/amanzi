/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Evaluators

*/

#ifndef AMANZI_EVALUATORS_POROSITY_COMPRESSIBLE_HH_
#define AMANZI_EVALUATORS_POROSITY_COMPRESSIBLE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Porosity.hh"

namespace Amanzi {
namespace Evaluators {

class Porosity_Compressible : public Porosity {
 public:
  explicit Porosity_Compressible(Teuchos::ParameterList& plist)
  {
    porosity_ = plist.get<double>("undeformed soil porosity");
    p_ref_ = plist.get<double>("reference pressure");
    c_ = plist.get<double>("pore compressibility");
    a_liquid_ = plist.get<double>("liquid thermal dilation", 0.0);
    a_rock_ = plist.get<double>("rock thermal dilation", 0.0);

    factor_ = porosity_ * c_;
  }
  ~Porosity_Compressible(){};

  // required methods from the base class
  virtual double PorosityValue(double p) override
  {
    double dp = p - p_ref_;
    return porosity_ * std::exp(c_ * dp);
  }
  virtual double dPorositydPressure(double p) override
  {
    double dp = p - p_ref_;
    return factor_ * std::exp(c_ * dp);
  }

  virtual double PorosityValueReference() override { return porosity_; }
  virtual std::pair<double, double> getThermalCoefficients() override
  {
    return std::make_pair(a_liquid_, a_rock_);
  }

 private:
  double porosity_, p_ref_, c_;
  double factor_;
  double a_liquid_, a_rock_;
};

} // namespace Evaluators
} // namespace Amanzi

#endif
