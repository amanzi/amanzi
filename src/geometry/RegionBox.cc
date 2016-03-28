/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A rectangular region in space, defined by two points

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
           Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
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
		     const Set_ID id,
                     const Point& p0,
		     const Point& p1,
                     const LifeCycleType lifecycle)
  : Region(name, id, true, BOX, p0.dim(), p0.dim(), lifecycle),
    p0_(p0),
    p1_(p1)
{
  if (p0_.dim() != p1_.dim()) {
    Errors::Message msg;
    msg << "Mismatch in dimensions of corner points of RegionBox \""
	<< Region::name() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  
  // Check if this is a reduced dimensionality box (e.g. even though
  // it is in 3D space it is a 2D box)

  int dim = p0.dim();
  for (int i=0; i!=p0.dim(); ++i)
    if (p0[i] == p1[i]) dim--;
  
  if (dim < p0.dim()) set_manifold_dimension(dim);
}

void
RegionBox::corners(Point *lo_corner, Point *hi_corner) const
{
  ASSERT(lo_corner != NULL);
  ASSERT(hi_corner != NULL);

  lo_corner->set(p0_);
  hi_corner->set(p1_);
}
  
// -------------------------------------------------------------
// RegionBox::between_
// -------------------------------------------------------------
bool
RegionBox::between_(const double& x, const double& x0, const double& x1)
{
  return  (x0+TOL >= x && x >= x1-TOL) || (x1+TOL >= x && x >= x0-TOL);
}

// -------------------------------------------------------------
// RegionBox::inside
// -------------------------------------------------------------
bool
RegionBox::inside(const Point& p) const
{

#ifdef ENABLE_DBC
  if (p.dim() != p0_.dim()) {
    Errors::Message msg;
    msg << "Mismatch in corner dimension of RegionBox \""
	<< name() << "\" and query point.";
    Exceptions::amanzi_throw(msg);
  }
#endif

  bool result(true);
  for (int i = 0; i!=p.dim(); ++i) {
    result = result && between_(p[i], p0_[i], p1_[i]);
  }
  return result;
}

// -------------------------------------------------------------
// RegionBox::is_degenerate (also indicate in how many dimensions)
// -------------------------------------------------------------
bool
RegionBox::is_degenerate(int *ndeg) const
{
  *ndeg = 0;
  for (int i=0; i!=p0_.dim(); ++i) {
    if (p0_[i] == p1_[i]) (*ndeg)++;    
  }

  if (*ndeg) return true;
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
