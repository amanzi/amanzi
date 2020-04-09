/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
RegionHalfSpace::RegionHalfSpace(const std::string& name, const int id,
                                 const Point& p, const Point& normal,
                                 const LifeCycleType lifecycle)
  : Region(name, id, true, PLANE, p.dim() - 1, p.dim(), lifecycle),
    p_(p),
    n_(normal / norm(normal))
{
  if (p_.dim() != n_.dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in point and normal dimensions of RegionHalfSpace "
         << Region::name();
    Exceptions::amanzi_throw(mesg);
  }
}


// -------------------------------------------------------------
// RegionHalfSpace::inside -- check if point is on plane
// -------------------------------------------------------------
bool
RegionHalfSpace::inside(const Point& p) const
{
#ifdef ENABLE_DBC
  if (p_.dim() != n_.dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in point dimension of RegionHalfSpace \""
         << Region::name() << "\" and query point.";
    Exceptions::amanzi_throw(mesg);
  }
#endif


  double res(0.0);

  for (int i = 0; i != p.dim(); ++i) { res += n_[i] * (p[i] - p_[i]); }

  return res > TOL ? false : true;
}

} // namespace AmanziGeometry
} // namespace Amanzi
