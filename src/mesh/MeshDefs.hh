// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
// file: MeshDefs.hh
// -------------------------------------------------------------
/**
 * @file   MeshDefs.hh
 * @author William A. Perkins
 * @date Mon May  2 13:03:23 2011
 * 
 * @brief  Various definitions needed by Mesh
 * 
 * 
 */
// -------------------------------------------------------------
// Created May  2, 2011 by William A. Perkins
// Last Change: Mon May  2 13:03:23 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _MeshDefs_hh_
#define _MeshDefs_hh_

#include <vector>

namespace Amanzi {
namespace AmanziMesh {

// Necessary typedefs and enumerations
      
typedef unsigned int Entity_ID;
typedef unsigned int Set_ID;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef std::vector<Set_ID> Set_ID_List;

// Cells (aka zones/elements) are the highest dimension entities in a mesh 
// Nodes (aka vertices) are lowest dimension entities in a mesh 
// Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
    
enum Entity_kind 
{
  NODE = 0,
  EDGE,
  FACE,
  CELL
};

// Check if Entity_kind is valid
inline 
bool entity_valid_kind (const Entity_kind kind) {
  return (kind >= NODE && kind <= CELL);
}

// Parallel status of entity 
    
enum Parallel_type 
{
  OWNED = 1,         // Owned by this processor
  GHOST = 2,         // Owned by another processor
  USED  = 3          // OWNED + GHOST
};

// Check if Parallel_type is valid

inline 
bool entity_valid_ptype (const Parallel_type ptype) {
  return (ptype >= OWNED && ptype <= USED);
}
    
// Standard element types and catchall (POLYGON/POLYHED)

enum Cell_type 
{
  UNKNOWN = 0,
  TRI = 1,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHED                // Polyhedron 
};
    
// Check if Cell_type is valid
inline 
bool cell_valid_type (const Cell_type type) {
  return (type >= TRI && type <= POLYHED); 
}

} // close namespace Amanzi 
} // close namespace AmanziMesh



#endif
