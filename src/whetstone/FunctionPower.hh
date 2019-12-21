/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Power function a (s - 0.5) ^ n.
*/

#ifndef AMANZI_WHETSTONE_FUNCTION_POWER_HH_
#define AMANZI_WHETSTONE_FUNCTION_POWER_HH_

#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class FunctionPower : public WhetStoneFunction {
 public:
  FunctionPower(double factor, int n) : factor_(factor), n_(n) {};
  ~FunctionPower() {};

  virtual double Value(const AmanziGeometry::Point& xp, double t) const override {
    return factor_ * std::pow(0.5 - t, n_);
  }

 private:
  double factor_;
  int n_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
