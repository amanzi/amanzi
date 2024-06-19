/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon ecoon@lanl.gov
*/

/* -------------------------------------------------------------------------
  Spline

 Cubic Hermite spline functor which is used for smoothing.  NOTE: this is NOT
 guaranteed to be monotonic!  In practice I hope it is usually!

------------------------------------------------------------------------- */

#include "Teuchos_SerialDenseMatrix.hpp"
#include "dbc.hh"

#include "Spline.hh"


namespace Amanzi {
namespace Utils {

void
Spline::Setup(double x1, double y1, double dy1, double x2, double y2, double dy2)
{
}


} // namespace Utils
} // namespace Amanzi
