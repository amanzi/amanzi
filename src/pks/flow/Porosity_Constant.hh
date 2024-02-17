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

#ifndef AMANZI_FLOW_POROSITY_CONSTANT_HH_
#define AMANZI_FLOW_POROSITY_CONSTANT_HH_

#include <utility>

#include "Teuchos_ParameterList.hpp"

#include "Porosity.hh"

namespace Amanzi {
namespace Flow {

class Porosity_Constant : public Porosity {
 public:
  explicit Porosity_Constant(Teuchos::ParameterList& plist)
  {
    porosity_ = plist.get<double>("value");
  }
  ~Porosity_Constant(){};

  // required methods from the base class
  virtual double PorosityValue(double p) override { return porosity_; }
  virtual double dPorositydPressure(double p) override { return 0.0; }
  virtual double getBiotCoefficient() override { return 1.0; }
  virtual std::pair<double, double> getThermalCoefficients() override
  {
    return std::make_pair(0.0, 0.0);
  }

 private:
  double porosity_;
};

} // namespace Flow
} // namespace Amanzi

#endif
