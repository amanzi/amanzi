/*
  A rectangular region in space, defined by two corner points and
  normals to its side.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

#include <vector>

#include "dbc.hh"
#include "errors.hh"

#include "Point.hh"
#include "RegionBoxVolumeFractions.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
RegionBoxVolumeFractions::RegionBoxVolumeFractions(
    const std::string& name, const Set_ID id,
    const Point& p0, const Point& p1,
    const std::vector<Point>& normals,
    const LifeCycleType lifecycle)
  : Region(name, id, true, BOX_VOF, p0.dim(), p0.dim(), lifecycle),
    p0_(p0),
    p1_(p1),
    normals_(normals)
{
  Errors::Message msg;
  if (p0_.dim() != p1_.dim()) {
    msg << "Mismatch in dimensions of corner points of RegionBoxVolumeFractions \""
	<< Region::name() << "\"";
    Exceptions::amanzi_throw(msg);

    for (int n = 0; n < normals.size(); ++n) {
      if (p0_.dim() != normals_[n].dim()) {
        msg << "Mismatch in dimensions of points and normals of RegionBoxVolumeFractions \""
	    << Region::name() << "\"";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  // Calculate the transformation tensor
  int dim = p0.dim();
  N_.set(dim);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      N_(j, i) = normals_[i][j];
    }
  }
  N_.Inverse();

  // Check if this is a reduced dimensionality box (e.g. even though
  // it is in 3D space it is a 2D box). Normalize the tensor to unit
  // reference cube.
  Point v1(dim);
  v1 = N_ * (p1_ - p0_);

  int mdim(dim);
  for (int i = 0; i != dim; ++i) {
    double len = v1[i];
    if (std::abs(len) < 1e-12 * norm(v1)) {
      mdim--;
    } else {
      for (int j = 0; j < dim; ++j) N_(i, j) /= len;
    }
  } 
  
  if (mdim < dim) set_manifold_dimension(mdim);
}


// -------------------------------------------------------------
// Implementation of virtual member function.
// -------------------------------------------------------------
bool RegionBoxVolumeFractions::inside(const Point& p) const
{
#ifdef ENABLE_DBC
  if (p.dim() != p0_.dim()) {
    Errors::Message msg;
    msg << "Mismatch in corner dimension of RegionBox \""
	<< name() << "\" and query point.";
    Exceptions::amanzi_throw(msg);
  }
#endif

  Point phat(N_ * (p - p0_));
 
  bool result(true);
  for (int i = 0; i != p.dim(); ++i) {
    result &= (phat[i] > -TOL && phat[i] < 1.0 + TOL);
  }
  return result;
}


// -------------------------------------------------------------
// Also indicates in how many dimensions the box is degenerate.
// -------------------------------------------------------------
bool RegionBoxVolumeFractions::is_degenerate(int *ndeg) const
{
  *ndeg = 0;
  for (int i = 0; i != p0_.dim(); ++i) {
    if (p0_[i] == p1_[i]) (*ndeg)++;    
  }

  if (*ndeg) return true;
  return false;
}

}  // namespace AmanziGeometry
}  // namespace Amanzi
