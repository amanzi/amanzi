/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Upwind function.
*/

#ifndef AMANZI_WHETSTONE_COMPOSITE_FUNCTION_HH_
#define AMANZI_WHETSTONE_COMPOSITE_FUNCTION_HH_

#include "Polynomial.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class FunctionComposite : public WhetStoneFunction {
 public:
  FunctionComposite(const WhetStoneFunction* f1, std::vector<const WhetStoneFunction*>& f2)
    : f1_(f1), f2_(f2){};
  ~FunctionComposite(){};

  virtual double Value(const AmanziGeometry::Point& xp) const
  {
    int d = f2_.size();
    AmanziGeometry::Point yp(d);
    for (int i = 0; i < d; ++i) yp[i] = f2_[i]->Value(xp);
    return f1_->Value(yp);
  }

 private:
  const WhetStoneFunction* f1_;
  std::vector<const WhetStoneFunction*> f2_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
