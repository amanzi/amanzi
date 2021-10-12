/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A planar (infinite) region in space, defined by a point and a normal.

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

#include "RegionPlane.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// RegionPlane:: constructors / destructor
// -------------------------------------------------------------
RegionPlane::RegionPlane(const std::string& name, 
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
    mesg << "Mismatch in point and normal dimensions of RegionPlane "
         << Region::name();
    Exceptions::amanzi_throw(mesg);
  }
}


// -------------------------------------------------------------
// RegionPlane::inside -- check if point is on plane
// -------------------------------------------------------------
bool
RegionPlane::inside(const Point& p) const
{
#ifdef ENABLE_DBC
  if (p.dim() != n_.dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in point dimension of RegionPlane \""
         << Region::name() << "\" and query point.";
    Exceptions::amanzi_throw(mesg);
  }
#endif

  double d(0.0), res(0.0);
  for (int i=0; i!=p.dim(); ++i) {
    res += n_[i]*p[i];
    d += n_[i]*p_[i];
  }
  res -= d;

  return std::abs(res) <= tol_ ? true : false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
