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

#include "SplinePolynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Calculate mesh velocity in cell c. 
****************************************************************** */
void SplinePolynomial::Setup(double x0, double f0, double df0,
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
};

}  // namespace WhetStone
}  // namespace Amanzi

