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
PointRegion::PointRegion(const Set_Name& name, 
			 const Set_ID id,
			 const Point& p,
                         const LifeCycleType lifecycle,
                         const VerboseObject *verbobj)
  : Region(name,id,0,lifecycle,verbobj), p_(p)
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


  if (p.dim() != p_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in dimension of PointRegion \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";

    const VerboseObject *verbobj = Region::verbosity_obj();
    if (verbobj && verbobj->os_OK(Teuchos::VERB_MEDIUM)) {
      Teuchos::OSTab tab = verbosity_obj()->getOSTab();
      *(verbobj->os()) << tempstr;
    }
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }

  bool result(true);

  for (int i = 0; i < p.dim(); ++i) 
    {
      result = result & (fabs(p[i]-p_[i]) < 1e-12);
    }

  return result;
}

} // namespace AmanziGeometry
} // namespace Amanzi
