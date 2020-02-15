/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Spline function is based on a one-dimensional polynomial which 
  provides easy derivatives and possiblity for generalization to
  multiple dimensions.
*/

#ifndef AMANZI_WHETSTONE_SPLINE_POLYNOMIAL_HH_
#define AMANZI_WHETSTONE_SPLINE_POLYNOMIAL_HH_

#include "Polynomial.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class SplinePolynomial : public WhetStoneFunction {
 public:
  SplinePolynomial(const AmanziGeometry::Point& x0, double f0, double df0,
                   const AmanziGeometry::Point& x1, double f1, double df1);
  ~SplinePolynomial() {};

  virtual double Value(const AmanziGeometry::Point& xp) const {
    return poly_.Value(xp);
  }

  // access
  const Polynomial poly() const { return poly_; } 

 private:
  Polynomial poly_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
