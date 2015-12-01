// Enumerated region is a user-provided list of mesh entities.
//
// Author: Ethan Coon
//

#include "EnumeratedSetRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {


// -------------------------------------------------------------
// EnumeratedSetRegion:: constructors / destructor
// -------------------------------------------------------------
EnumeratedSetRegion::EnumeratedSetRegion(const Set_Name& name, 
				   const Set_ID id,
				   const std::string& entity_str,
				   const Entity_ID_List& ents,
                                   const LifeCycleType lifecycle,
                                   const VerboseObject *verbobj)
  : Region(name,id,3,lifecycle,verbobj),
    entity_str_(entity_str),
    entities_(ents)    
{
  // empty
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension
}

EnumeratedSetRegion::EnumeratedSetRegion(const EnumeratedSetRegion& old)
  : Region(old),
    entity_str_(old.entity_str_),
    entities_(old.entities_)
{}

EnumeratedSetRegion::~EnumeratedSetRegion(void)
{}


// -------------------------------------------------------------
// EnumeratedSetRegion::inside
// -------------------------------------------------------------
bool
EnumeratedSetRegion::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for enumerated sets");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
