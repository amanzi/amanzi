/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
      Rao Garimella
      Ethan Coon (ecoon@lanl.gov)
*/

/*
  A rectangular region in space, defined by two opposite corners.

*/

#include "dbc.hh"
#include "errors.hh"

#include "RegionBox.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
RegionBox::RegionBox(const std::string& name,
                     const int id,
                     const Point& p0,
                     const Point& p1,
                     const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::BOX, p0.dim(), p0.dim(), lifecycle), p0_(p0), p1_(p1)
{
  if (p0_.dim() != p1_.dim()) {
    Errors::Message msg;
    msg << "Mismatch in dimensions of corner points of RegionBox \"" << Region::get_name() << "\"";
    Exceptions::amanzi_throw(msg);
  }

  // Check if this is a reduced dimensionality box (e.g. even though
  // it is in 3D space it is a 2D box)
  int dim = p0.dim();
  for (int i = 0; i != p0.dim(); ++i)
    if (p0[i] == p1[i]) dim--;

  if (dim < p0.dim()) set_manifold_dimension(dim);

  // create opposite corners
  for (int i = 0; i != p0.dim(); ++i) {
    p0_[i] = std::min(p0[i], p1[i]);
    p1_[i] = std::max(p0[i], p1[i]);
  }
}


// -------------------------------------------------------------
// RegionBox::inside
// -------------------------------------------------------------
bool
RegionBox::inside(const Point& p) const
{
  // #ifdef ENABLE_DBC
  if (p.dim() != p0_.dim()) {
    Errors::Message msg;
    msg << "Mismatch in corner dimension of RegionBox \"" << get_name() << "\" and query point.";
    Exceptions::amanzi_throw(msg);
  }
  // #endif

  for (int i = 0; i != p.dim(); ++i) {
    if (p[i] < p0_[i] - TOL) return false;
    if (p[i] > p1_[i] + TOL) return false;
  }
  return true;
}


// -------------------------------------------------------------
// RegionBox::is_degenerate (also indicate in how many dimensions)
// -------------------------------------------------------------
bool
RegionBox::is_degenerate(int* ndeg) const
{
  *ndeg = 0;
  for (int i = 0; i != p0_.dim(); ++i) {
    if (p0_[i] == p1_[i]) (*ndeg)++;
  }

  if (*ndeg) return true;
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
