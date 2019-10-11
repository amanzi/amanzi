/*
  WhetStone, Version 2.2
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
  FunctionUpwindPlus(const WhetStoneFunction* f) : f_(f){};
  ~FunctionUpwindPlus(){};

  virtual double Value(const AmanziGeometry::Point& xp) const
  {
    return std::max(f_->Value(xp), 0.0);
  }

 private:
  const WhetStoneFunction* f_;
};


class FunctionUpwindMinus : public WhetStoneFunction {
 public:
  FunctionUpwindMinus(const WhetStoneFunction* f) : f_(f){};
  ~FunctionUpwindMinus(){};

  virtual double Value(const AmanziGeometry::Point& xp) const
  {
    return std::min(f_->Value(xp), 0.0);
  }

 private:
  const WhetStoneFunction* f_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
