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

namespace Amanzi {
namespace WhetStone {

// base class for splines
class SplinePolynomial {
 public:
  SplinePolynomial(){};
  ~SplinePolynomial(){};

  virtual void Setup(double x0, double f0, double df0, double x1, double f1, double df1) = 0;

  virtual double Value(double x) const = 0;
  virtual double GradientValue(double x) const = 0;
};


// cubic interpolant between two points
class SplineCubic : public SplinePolynomial {
 public:
  SplineCubic(){};
  ~SplineCubic(){};

  virtual void Setup(double x0, double f0, double df0, double x1, double f1, double df1);

  virtual double Value(double x) const;
  virtual double GradientValue(double x) const;

  // access
  const Polynomial poly() const { return poly_; }
  const Polynomial grad() const { return grad_; }

 private:
  Polynomial poly_;
  Polynomial grad_;
};


// linear interpolant exterior to the interval defined by two points
class SplineExteriorLinear : public SplinePolynomial {
 public:
  SplineExteriorLinear(){};
  ~SplineExteriorLinear(){};

  virtual void Setup(double x0, double f0, double df0, double x1, double f1, double df1);

  virtual double Value(double x) const;
  virtual double GradientValue(double x) const;

 private:
  double x0_, f0_, df0_, x1_, f1_, df1_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
