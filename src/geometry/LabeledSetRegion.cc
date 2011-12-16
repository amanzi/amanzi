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

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class LabeledSetRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// LabeledSetRegion:: constructors / destructor
// -------------------------------------------------------------
LabeledSetRegion::LabeledSetRegion(const std::string name, 
				   const unsigned int id,
				   const std::string entity_str,
                                   const std::string file,
                                   const std::string format,
                                   const std::string label)
  : Region(name,id,3),entity_str_(entity_str),
    file_(file), format_(format), label_(label)
{
  // empty
  // Region dimension is set arbitrarily as 3 since the set of
  // entities in the mesh will determine the dimension
}

LabeledSetRegion::LabeledSetRegion(const char *name, 
                                   const unsigned int id,
				   const std::string entity_str,
                                   const std::string file,
                                   const std::string format,
                                   const std::string label)
  : Region(name,id,3),entity_str_(entity_str),
    file_(file), format_(format), label_(label)
{
  // empty
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
  std::cerr << "In/out check not implemented for labeled sets" << std::endl;
  throw std::exception();
}

} // namespace AmanziGeometry
} // namespace Amanzi
