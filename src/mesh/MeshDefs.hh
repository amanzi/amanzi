/*
  Mesh

  Copyright 2010-201x held jointly by LANL, ORNL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William A. Perkins
           Rao Garimella

 Various definitions needed by Mesh class.
*/

#ifndef AMANZI_MESH_DEFS_HH_
#define AMANZI_MESH_DEFS_HH_

#include <vector>
#include <string>
#include "boost/algorithm/string.hpp"

#include "errors.hh"
#include "Point.hh"

namespace Amanzi {
namespace AmanziMesh {

// Necessary typedefs and enumerations
typedef int Entity_ID;
typedef int Entity_GID;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef std::vector<Entity_GID> Entity_GID_List;
typedef std::vector<int> Entity_Direction_List;
typedef std::vector<AmanziGeometry::Point> Point_List;
typedef int Set_ID;

// Cells (aka zones/elements) are the highest dimension entities in a mesh
// Nodes (aka vertices) are lowest dimension entities in a mesh
// Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
// Entity_kind::BOUNDARY_FACE is a special type of entity that is need so that process
// kernels can define composite vectors (see src/data_structures) on
// exterior boundary faces of the mesh only
enum class Entity_kind
{
  UNKNOWN = 0,
  NODE = 1,
  EDGE = 2,
  FACE = 3,
  CELL = 4,
  BOUNDARY_FACE,
};

// Check if Entity_kind is valid
// inline
// bool validEntityKind (const int kind) {
//   return (kind >= Entity_kind::NODE && kind <= Entity_kind::CELL);
// }

// entity kind from string
inline
Entity_kind createEntityKind(const std::string& instring)
{
  std::string estring = instring; // note not done in signature to throw a better error
  boost::algorithm::to_lower(estring);
  if (estring == "cell") return Entity_kind::CELL;
  else if (estring == "face") return Entity_kind::FACE;
  else if (estring == "boundary_face") return Entity_kind::BOUNDARY_FACE;
  else if (estring == "edge") return Entity_kind::EDGE;
  else if (estring == "node") return Entity_kind::NODE;
  else {
    Errors::Message msg;
    msg << "Unknown entity kind string: \"" << instring << "\", valid are \"cell\", \"face\", \"boundary_face\", \"edge\", and \"node\".";
    Exceptions::amanzi_throw(msg);
    return Entity_kind::NODE;
  }
}

// string from entity kind
inline
std::string to_string(const Entity_kind kind)
{
  switch(kind) {
    case(Entity_kind::CELL): return "cell";
    case(Entity_kind::FACE): return "face";
    case(Entity_kind::BOUNDARY_FACE): return "boundary_face";
    case(Entity_kind::EDGE): return "edge";
    case(Entity_kind::NODE): return "node";
    default: return "unknown";
  }
}

// Parallel status of entity
enum class Parallel_type {
  UNKNOWN = 0,
  OWNED = 1,  // Owned by this processor
  GHOST = 2,  // Owned by another processor
  ALL = 3     // OWNED + GHOST
};

// Standard element types and catchall (POLYGON/POLYHED)
enum class Cell_type {
  UNKNOWN = 0,
  TRI = 1,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHED  // Polyhedron
};

// // Check if Cell_type is valid
// inline
// bool validCellType (const int cell_type) {
//   return (cell_type >= Cell_type::TRI && cell_type <= Cell_type::POLYHED);
// }

// Types of partitioners (partitioning scheme bundled into the name)
enum class Partitioner_type {
  METIS = 0, // default
  ZOLTAN_GRAPH,
  ZOLTAN_RCB
};

// Return an string description for each partitioner type
inline
std::string to_string(const Partitioner_type partitioner_type) {
  switch(partitioner_type) {
    case(Partitioner_type::METIS): return "Partitioner_type::METIS";
    case(Partitioner_type::ZOLTAN_GRAPH): return "Partitioner_type::ZOLTAN_GRAPH";
    case(Partitioner_type::ZOLTAN_RCB): return "Partitioner_type::ZOLTAN_RCB";
    default: return "unknown";
  }
}

inline
Partitioner_type createPartitionerType(const std::string& pstring) {
  if (pstring == "metis" || pstring == "METIS") {
    return Partitioner_type::METIS;
  } else if (pstring == "ZOLTAN_GRAPH" || pstring == "zoltan_graph") {
    return Partitioner_type::ZOLTAN_GRAPH;
  } else if (pstring == "ZOLTAN_RCB" || pstring == "zoltan_rcb") {
    return Partitioner_type::ZOLTAN_RCB;
  } else {
    Errors::Message msg;
    msg << "Unknown Partitioner_type string: \"" << pstring << "\", valid are \"metis\", \"zoltan_graph\", \"zoltan_rcb\"";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziMesh
}  // namespace Amanzi

#endif
