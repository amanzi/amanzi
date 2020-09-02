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

#include "errors.hh"

#include "SplinePolynomial.hh"
#include "VectorObjects.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Setup cubic polynomial
****************************************************************** */
void SplineCubic::Setup(double x0, double f0, double df0,
                        double x1, double f1, double df1)
{
  double dx = x1 - x0;
  double df = f1 - f0;

  AmanziGeometry::Point p0(1);
  p0[0] = x0;

  poly_.Reshape(1, 3);
  poly_.set_origin(p0);

  poly_(0) = f0;
  poly_(1) = df0;
  poly_(2) = (3 * df - dx * (2 * df0 + df1)) / (dx * dx);
  poly_(3) = (-2 * df + dx * (df0 + df1)) / (dx * dx * dx);

  grad_ = (Gradient(poly_))[0];
};


/* ******************************************************************
* Calculate value
****************************************************************** */
double SplineCubic::Value(double x) const
{
  AmanziGeometry::Point xp(1);
  xp[0] = x;
  return poly_.Value(xp);
}


/* ******************************************************************
* Calculate value of the gradient
****************************************************************** */
double SplineCubic::GradientValue(double x) const
{
  AmanziGeometry::Point xp(1);
  xp[0] = x;
  return grad_.Value(xp);
}


/* ******************************************************************
* Setup of linear exterior interpolants
****************************************************************** */
void SplineExteriorLinear::Setup(double x0, double f0, double df0,
                                 double x1, double f1, double df1)
{
  x0_ = x0;
  f0_ = f0;
  df0_ = df0;

  x1_ = x1;
  f1_ = f1;
  df1_ = df1;
};


/* ******************************************************************
* Calculate value
****************************************************************** */
double SplineExteriorLinear::Value(double x) const
{
  if (x <= x0_) return f0_ + df0_ * (x - x0_);
  if (x >= x1_) return f1_ + df1_ * (x - x1_);

  Errors::Message msg("Evaluation of exterior interpolant outside its domain.");
  Exceptions::amanzi_throw(msg);
  return 0.0;
}


/* ******************************************************************
* Calculate gradient value
****************************************************************** */
double SplineExteriorLinear::GradientValue(double x) const
{
  if (x <= x0_) return df0_;
  if (x >= x1_) return df1_;

  Errors::Message msg("Evaluation of exterior interpolant outside its domain.");
  Exceptions::amanzi_throw(msg);
  return 0.0;
}

}  // namespace WhetStone
}  // namespace Amanzi

