/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/**
 * @file   LabeledSetRegion.hh
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Declaration of Labeled Set Region class which derives 
 *         its definition from a named set of mesh entities in
 *         a mesh file
 * 
 * 
 */

#ifndef _LabeledSetRegion_hh_
#define _LabeledSetRegion_hh_

#include "Region.hh"

namespace Amanzi {

  namespace AmanziGeometry {

// -------------------------------------------------------------
//  class LabeledSetRegion
// -------------------------------------------------------------
/// A region defined by a set of mesh entities in a mesh file
///
/// Strictly speaking, we should tie this region class to a particular
/// mesh or mesh file but that cause a circular dependency of meshes
/// on regions and of labeled set regions on meshes. We will rely on the 
/// fact that when a mesh is created specifying a geometric model, it
/// will create mesh entity sets based on the labeled sets in that 
/// geometric model. 
///
/// If we need to change this behavior, then we can make a forward
/// declaration of AmanziMesh::Mesh, make the Mesh class a friend, add
/// a mesh variable to this class and have a protected method to set
/// the mesh

class LabeledSetRegion : public Region {
public:

  /// Default constructor 

  LabeledSetRegion(const std::string name, 
                   const unsigned int id, 
                   const std::string entity_str,
                   const std::string file,
                   const std::string format,
                   const std::string label);

  LabeledSetRegion(const char *name, 
                   const unsigned int id, 
                   const std::string entity_str,
                   const std::string file,
                   const std::string format,
                   const std::string label);


  /// Protected copy constructor to avoid unwanted copies.
  LabeledSetRegion(const LabeledSetRegion& old);

  /// Destructor
  ~LabeledSetRegion(void);

  // Type of the region
  inline RegionType type() const { return LABELEDSET; }

  // Label in the file
  inline std::string label() const { return label_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  inline std::string entity_str() const { return entity_str_; }

protected:  
  const std::string entity_str_; // what kind of entities make up this set
  const std::string file_; // which file are we supposed to read it from
  const std::string format_; // format of the file
  const std::string label_; // Label used to identify set in the file (may be different from name)
};

/// A smart pointer to LabeledSetRegion instances
// typedef Teuchos::RCP<LabeledSetRegion> LabeledSetRegionPtr;

// RVG: I am not able to correctly code a region factory using smart
// pointers so I will revert to a simpler definition

typedef LabeledSetRegion *LabeledSetRegionPtr;

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
