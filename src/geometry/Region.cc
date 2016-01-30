/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   Region.cc
 * @author William A. Perkins
 * @date Mon Aug  1 10:05:25 2011
 * 
 * @brief  
 * 
 * 
 */

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class Region
// -------------------------------------------------------------

// -------------------------------------------------------------
// Region:: constructors / destructor
// -------------------------------------------------------------
Region::Region() :   verbosity_obj_(NULL)
{
  name_ = "";
  id_ = 0;
}

Region::Region(const Region& old) :
  topo_dimension_(old.dimension()), name_(old.name()), id_(old.id()),
  lifecycle_(old.lifecycle()), verbosity_obj_(old.verbosity_obj())
{
  // empty
}


Region::Region(const Set_Name& name, const Set_ID id, 
               const unsigned int dim, const LifeCycleType lifecycle,
               const VerboseObject *verbobj) :
  name_(name), id_(id), topo_dimension_(dim), lifecycle_(lifecycle),
  verbosity_obj_(verbobj)
{
}

Region::~Region(void)
{
  // empty
}

// Get the extents of the Region

// void Region::extents(Point *pmin, Point *pmax) const
// {
//   pmin = new Point(min_pnt);
//   pmax = new Point(max_pnt);
// }

} // namespace AmanziGeometry
} // namespace Amanzi
