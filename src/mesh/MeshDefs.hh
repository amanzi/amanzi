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

#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"

#include "errors.hh"
#include "AmanziTypes.hh"
#include "Point.hh"

#include "MeshView.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// Typedefs
//
using size_type = Kokkos::MeshView<int*, Kokkos::DefaultHostExecutionSpace>::size_type;
using Entity_ID = int;
using Entity_GID = int;
using Set_ID = int;
using Direction_type = int;

namespace Impl {
template<MemSpace_kind MEM>
struct MemorySpace {
  using space = Amanzi::DefaultMemorySpace;
};

template<>
struct MemorySpace<MemSpace_kind::HOST> {
  using space = Amanzi::DefaultHostMemorySpace;
};

} // namespace Impl

//
// Views are on host or device
//
template<typename T, MemSpace_kind MEM=MemSpace_kind::HOST>
using View_type = Kokkos::MeshView<T*, typename Impl::MemorySpace<MEM>::space>;

using Entity_ID_View = View_type<Entity_ID>;
using cEntity_ID_View = View_type<const Entity_ID>;
using Entity_GID_View = View_type<Entity_GID>;
using cEntity_GID_View = View_type<const Entity_GID>;
using Entity_Direction_View = View_type<int>;
using cEntity_Direction_View = View_type<const int>;
using Point_View = View_type<AmanziGeometry::Point>;
using cPoint_View = View_type<const AmanziGeometry::Point>;
using Double_View = View_type<double>;
using cDouble_View = View_type<const double>;

//
// View for host only
//
using Entity_ID_List = std::vector<Entity_ID>;
using Entity_GID_List = std::vector<Entity_ID>;
using Entity_Direction_List = std::vector<int>;
using Point_List = std::vector<AmanziGeometry::Point>;
using Double_List = std::vector<double>;
template<typename T> using RaggedArray_List = std::vector<std::vector<T>>;


template<typename T> using DualView_type = Kokkos::MeshDualView<T*, Kokkos::DefaultHostExecutionSpace>;
using Entity_ID_DualView = DualView_type<Entity_ID>;
using Entity_GID_DualView = DualView_type<Entity_GID>;
using Entity_Direction_DualView = DualView_type<int>;
using Point_DualView = DualView_type<AmanziGeometry::Point>;
using Double_DualView = DualView_type<double>;

//
// This may not be necessary?
//
struct RaggedArray_View {
  Entity_ID_View rows;
  Entity_ID_View entries;
};


// Cells (aka zones/elements) are the highest dimension entities in a mesh
// Nodes (aka vertices) are lowest dimension entities in a mesh
// Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
// Entity_kind::BOUNDARY_FACE is a special type of entity that is need so that process
// kernels can define composite vectors (see src/data_structures) on
// exterior boundary faces of the mesh only
enum Entity_kind : int
{
  UNKNOWN = 0,
  NODE = 1,
  EDGE = 2,
  FACE = 3,
  CELL = 4,
  BOUNDARY_NODE = 11,
  BOUNDARY_FACE = 13
};


// entity kind from string
inline
Entity_kind createEntityKind(const std::string& instring)
{
  std::string estring = instring; // note not done in signature to throw a better error
  transform(estring.begin(), estring.end(), estring.begin(), ::tolower);
  if (estring == "cell") return Entity_kind::CELL;
  else if (estring == "face") return Entity_kind::FACE;
  else if (estring == "boundary_face") return Entity_kind::BOUNDARY_FACE;
  else if (estring == "boundary_node") return Entity_kind::BOUNDARY_NODE;
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
    case(Entity_kind::BOUNDARY_NODE): return "boundary_node";
    case(Entity_kind::EDGE): return "edge";
    case(Entity_kind::NODE): return "node";
    default: return "unknown";
  }
}

// Parallel status of entity
enum class Parallel_kind {
  UNKNOWN = 0,
  OWNED = 1,  // Owned by this processor
  GHOST = 2,  // Owned by another processor
  ALL = 3     // OWNED + GHOST
};

inline
std::string to_string(const Parallel_kind ptype)
{
  switch(ptype) {
    case(Parallel_kind::UNKNOWN): return "UNKNOWN";
    case(Parallel_kind::OWNED): return "OWNED";
    case(Parallel_kind::GHOST): return "GHOST";
    case(Parallel_kind::ALL): return "ALL";
    default: return "unknown";
  }
}

// Standard element types and catchall (POLYGON/POLYHED)
enum class Cell_kind {
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

// string from entity kind
inline
std::string to_string(const Cell_kind ctype)
{
  switch(ctype) {
    case(Cell_kind::TRI): return "cell type: triangle";
    case(Cell_kind::QUAD): return "cell type: quadrilateral";
    case(Cell_kind::POLYGON): return "cell type: polygon";
    case(Cell_kind::TET): return "cell type: tetrahedron";
    case(Cell_kind::PRISM): return "cell type: prism";
    case(Cell_kind::PYRAMID): return "cell type: pyramid";
    case(Cell_kind::HEX): return "cell type: hexahedron";
    case(Cell_kind::POLYHED): return "cell type: polyhedron";
    default: return "cell type: unknown";
  }
}


// Types of partitioners (partitioning scheme bundled into the name)
enum class Partitioner_kind {
  METIS = 0, // default
  ZOLTAN_GRAPH,
  ZOLTAN_RCB
};

// Return an string description for each partitioner type
inline
std::string to_string(const Partitioner_kind partitioner_type) {
  switch(partitioner_type) {
    case(Partitioner_kind::METIS): return "Partitioner_kind::METIS";
    case(Partitioner_kind::ZOLTAN_GRAPH): return "Partitioner_kind::ZOLTAN_GRAPH";
    case(Partitioner_kind::ZOLTAN_RCB): return "Partitioner_kind::ZOLTAN_RCB";
    default: return "unknown";
  }
}

inline
Partitioner_kind createPartitionerType(const std::string& pstring) {
  if (pstring == "metis" || pstring == "METIS") {
    return Partitioner_kind::METIS;
  } else if (pstring == "ZOLTAN_GRAPH" || pstring == "zoltan_graph") {
    return Partitioner_kind::ZOLTAN_GRAPH;
  } else if (pstring == "ZOLTAN_RCB" || pstring == "zoltan_rcb") {
    return Partitioner_kind::ZOLTAN_RCB;
  } else {
    Errors::Message msg;
    msg << "Unknown Partitioner_kind string: \"" << pstring << "\", valid are \"metis\", \"zoltan_graph\", \"zoltan_rcb\"";
    Exceptions::amanzi_throw(msg);
  }
  return Partitioner_kind::METIS;
}

enum class AccessPattern_kind {
  DEFAULT,
  ANY,
  CACHE,
  COMPUTE,
  FRAMEWORK
};

using MeshSets = std::map<std::tuple<std::string,Entity_kind,Parallel_kind>, Entity_ID_DualView>;
using MeshSetVolumeFractions = std::map<std::tuple<std::string,Entity_kind,Parallel_kind>, Double_DualView>;


}  // namespace AmanziMesh
}  // namespace Amanzi


