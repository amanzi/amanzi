/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William A. Perkins
      Rao Garimella
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_MESH_DEFS_HH_
#define AMANZI_MESH_DEFS_HH_

#include <vector>
#include <string>
#include "boost/algorithm/string.hpp"

#include "AmanziTypes.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {

// Necessary typedefs and enumerations
typedef int Set_ID;
typedef LO Entity_ID;
typedef std::vector<Entity_ID> Entity_ID_List;
typedef Kokkos::View<Entity_ID*> Entity_ID_View;

// Recongnize special meshes
enum Mesh_type {
  RECTANGULAR, // Equivalent of structured but can't use i,j,k notation
  GENERAL      // general unstructured
};

// Cells (aka zones/elements) are the highest dimension entities in a mesh
// Nodes (aka vertices) are lowest dimension entities in a mesh
// Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
// BOUNDARY_FACE is a special type of entity that is need so that process
// kernels can define composite vectors (see src/data_structures) on
// exterior boundary faces of the mesh only
enum Entity_kind { NODE = 0, EDGE, FACE, CELL, BOUNDARY_FACE, UNKNOWN };

// Check if Entity_kind is valid
inline bool
entity_valid_kind(const Entity_kind kind)
{
  return (kind >= NODE && kind <= CELL);
}

// entity kind from string
inline Entity_kind
entity_kind(const std::string& instring)
{
  std::string estring =
    instring; // note not done in signature to throw a better error
  boost::algorithm::to_lower(estring);
  if (estring == "cell")
    return CELL;
  else if (estring == "face")
    return FACE;
  else if (estring == "boundary_face")
    return BOUNDARY_FACE;
  else if (estring == "edge")
    return EDGE;
  else if (estring == "node")
    return NODE;
  else {
    Errors::Message msg;
    msg << "Unknown entity kind string: \"" << instring
        << "\", valid are \"cell\", \"face\", \"boundary_face\", \"edge\", and "
           "\"node\".";
    Exceptions::amanzi_throw(msg);
    return NODE;
  }
}

// string from entity kind
inline std::string
entity_kind_string(Entity_kind kind)
{
  switch (kind) {
  case (CELL):
    return "cell";
  case (FACE):
    return "face";
  case (BOUNDARY_FACE):
    return "boundary_face";
  case (EDGE):
    return "edge";
  case (NODE):
    return "node";
  default:
    return "unknown";
  }
}

// Parallel status of entity
enum class Parallel_type : int {
  PTYPE_UNKNOWN = 0,
  OWNED = 1,             // Owned by this processor
  GHOST = 2,             // Owned by another processor
  ALL = 3,               // OWNED + GHOST
  PARALLEL_TYPE_SIZE = 4 // Keep this element last for the size of enum
};

// Check if Parallel_type is valid
inline bool
entity_valid_ptype(const Parallel_type ptype)
{
  return (ptype >= Parallel_type::OWNED && ptype <= Parallel_type::ALL);
}

// Standard element types and catchall (POLYGON/POLYHED)
enum Cell_type {
  CELLTYPE_UNKNOWN = 0,
  TRI = 1,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHED // Polyhedron
};

// Check if Cell_type is valid
inline bool
cell_valid_type(const Cell_type type)
{
  return (type >= TRI && type <= POLYHED);
}

// Types of partitioners (partitioning scheme bundled into the name)
enum class Partitioner_type : std::uint8_t { METIS, ZOLTAN_GRAPH, ZOLTAN_RCB };

constexpr int NUM_PARTITIONER_TYPES = 3;
constexpr Partitioner_type PARTITIONER_DEFAULT = Partitioner_type::METIS;

// Return an string description for each partitioner type
inline std::string
Partitioner_type_string(const Partitioner_type partitioner_type)
{
  static std::string partitioner_type_str[NUM_PARTITIONER_TYPES] = {
    "Partitioner_type::METIS",
    "Partitioner_type::ZOLTAN_GRAPH",
    "Partitioner_type::ZOLTAN_RCB"
  };

  int iptype = static_cast<int>(partitioner_type);
  return (iptype >= 0 && iptype < NUM_PARTITIONER_TYPES) ?
           partitioner_type_str[iptype] :
           "";
}

// Output operator for Partitioner_type
inline std::ostream&
operator<<(std::ostream& os, const Partitioner_type& partitioner_type)
{
  os << " " << Partitioner_type_string(partitioner_type) << " ";
  return os;
}

// Types of partitioning algorithms - Add as needed in the format METIS_RCB etc.
enum class Partitioning_scheme { DEFAULT };

} // namespace AmanziMesh
} // namespace Amanzi

#endif
