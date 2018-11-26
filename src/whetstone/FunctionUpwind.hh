/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Upwind function.
*/

#ifndef AMANZI_WHETSTONE_UPWIND_FUNCTION_HH_
#define AMANZI_WHETSTONE_UPWIND_FUNCTION_HH_

#include "Polynomial.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class FunctionUpwindPlus : public WhetStoneFunction {
 public:
  FunctionUpwindPlus(const WhetStoneFunction* un, const WhetStoneFunction* f) : un_(un), f_(f) {};
  ~FunctionUpwindPlus() {};

  virtual double Value(const AmanziGeometry::Point& xp) const {
    return (un_->Value(xp) >= 0.0) * f_->Value(xp);
  }

 private:
  const WhetStoneFunction *un_, *f_;
};


class FunctionUpwindMinus : public WhetStoneFunction {
 public:
  FunctionUpwindMinus(const WhetStoneFunction* un, const WhetStoneFunction* f) : un_(un), f_(f) {};
  ~FunctionUpwindMinus() {};

  virtual double Value(const AmanziGeometry::Point& xp) const {
    return (un_->Value(xp) <= 0.0) * f_->Value(xp);
  }

 private:
  const WhetStoneFunction *un_, *f_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
