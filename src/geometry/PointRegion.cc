/**
 * @file   PointRegion.cc
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Implementation of PointRegion class
 * 
 * 
 */

#include "PointRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class PointRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// PointRegion:: constructors / destructor
// -------------------------------------------------------------
PointRegion::PointRegion(const std::string name, 
			 const unsigned int id,
			 const Point& p)
  : Region(name,id,0), p_(p)
{
}

PointRegion::PointRegion(const char *name, const unsigned int id,
			 const Point& p)
  : Region(name,id,0), p_(p)
{
}

PointRegion::PointRegion(const PointRegion& old)
  : Region(old), p_(old.p_)
{
  // empty
}

PointRegion::~PointRegion(void)
{
  
}

// -------------------------------------------------------------
// PointRegion::inside -- check if input point is coincident with point
// -------------------------------------------------------------
bool
PointRegion::inside(const Point& p) const
{

#ifdef ENABLE_DBC
  if (p.dim() != p_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in dimension of PointRegion \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif

  bool result(true);

  for (int i = 0; i < p.dim(); ++i) 
    {
      result = result & (fabs(p[i]-p_[i]) < 1e-12);
    }

  return result;
}

} // namespace AmanziGeometry
} // namespace Amanzi
