/*
  A halfspace (infinite) region in space, defined by a plane.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
  Ethan Coon (ecoon@lanl.gov)
*/

#include <cmath>

#include "dbc.hh"
#include "errors.hh"
#include "Point.hh"

#include "RegionHalfSpace.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// RegionHalfSpace:: constructors / destructor
// -------------------------------------------------------------
RegionHalfSpace::RegionHalfSpace(const std::string& name, 
                                 const int id,
                                 const Point& p,
                                 const Point& normal,
                                 const LifeCycleType lifecycle)
    : Region(name, id, true, PLANE, p.dim()-1, p.dim(), lifecycle),
      p_(p),
      n_(normal/norm(normal))
{
  if (p_.dim() != n_.dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in point and normal dimensions of RegionHalfSpace "
         << Region::get_name();
    Exceptions::amanzi_throw(mesg);
  }
}


// -------------------------------------------------------------
// RegionHalfSpace::inside -- check if point is on plane
// -------------------------------------------------------------
bool
RegionHalfSpace::inside(const Point& p) const
{
  double res(0.0);
 
  for (int i = 0; i != p.dim(); ++i) {
    res += n_[i] * (p[i] - p_[i]);
  }

  return res > tol_ ? false : true;
}

} // namespace AmanziGeometry
} // namespace Amanzi
