/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

// Enumerated region is a user-provided list of mesh entities.
//
// I'm not particularly happy with this design -- it would be nice to
// have a way of saying, "this set is only valid on this mesh".  Maybe
// this is true of all regions, however. --etc
//
// Author: Ethan Coon
//

#ifndef AMANZI_ENUMERATED_REGION_HH_
#define AMANZI_ENUMERATED_REGION_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class EnumeratedSetRegion
// -------------------------------------------------------------
/// A rectangular region in space, defined by two points

class EnumeratedSetRegion : public Region {
public:

  /// Default constructor uses two corner points (order not important).

  EnumeratedSetRegion(const Set_Name& name,
                      const Set_ID id,
                      const std::string& entity_str,
                      const Entity_ID_List& ents,
                      const LifeCycleType lifecycle=PERMANENT,
                      const VerboseObject *verbobj=NULL);

  /// Protected copy constructor to avoid unwanted copies.
  EnumeratedSetRegion(const EnumeratedSetRegion& old);

  /// Destructor
  ~EnumeratedSetRegion(void);


  // Type of the region
  inline
  RegionType type() const { return ENUMERATEDSET; }

  /// Get the first point defining the region
  inline
  const std::vector<Entity_ID>& entities() const { return entities_; }

  inline std::string entity_str() const { return entity_str_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:
  const std::string entity_str_; // what kind of entities make up this set
  const Entity_ID_List entities_; // list of those included
};

/// A smart pointer to EnumeratedSetRegion instances
//
// typedef Teuchos::RCP<EnumeratedSetRegion> EnumeratedSetRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef EnumeratedSetRegion *EnumeratedSetRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
