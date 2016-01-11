/**
 * @file   LabeledSetRegion.cc
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Implementation of Labeled Set Region class which derives its
 *         definition from named set of mesh entities in a mesh file
 * 
 * 
 */

#include "LabeledSetRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class LabeledSetRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// LabeledSetRegion:: constructors / destructor
// -------------------------------------------------------------
LabeledSetRegion::LabeledSetRegion(const Set_Name& name, 
				   const Set_ID id,
				   const std::string& entity_str,
                                   const std::string& file,
                                   const std::string& format,
                                   const std::string& label,
                                   const LifeCycleType lifecycle,
                                   const VerboseObject *verbobj)
  : Region(name,id,3,lifecycle,verbobj),entity_str_(entity_str),
    file_(file), format_(format), label_(label)
{
  // empty
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension
}

LabeledSetRegion::LabeledSetRegion(const LabeledSetRegion& old)
  : Region(old)
{
  // empty
}

LabeledSetRegion::~LabeledSetRegion(void)
{
  // empty
}


// -------------------------------------------------------------
// LabeledSetRegion::inside
// -------------------------------------------------------------
bool
LabeledSetRegion::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for labeled sets");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
