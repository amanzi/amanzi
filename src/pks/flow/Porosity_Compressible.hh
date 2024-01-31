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

#ifndef AMANZI_FLOW_POROSITY_COMPRESSIBLE_HH_
#define AMANZI_FLOW_POROSITY_COMPRESSIBLE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Porosity.hh"

namespace Amanzi {
namespace Flow {

class Porosity_Compressible : public Porosity {
 public:
  explicit Porosity_Compressible(Teuchos::ParameterList& plist)
  {
    porosity_ = plist.get<double>("undeformed soil porosity");
    p_ref_ = plist.get<double>("reference pressure");
    c_ = plist.get<double>("pore compressibility");
    b_ = plist.get<double>("biot coefficient", 1.0);

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

  virtual double getBiotCoefficient() override { return b_; }

 private:
  double porosity_, p_ref_, c_;
  double factor_;
  double b_;
};

} // namespace Flow
} // namespace Amanzi

#endif
