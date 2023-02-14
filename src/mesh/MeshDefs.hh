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
#include "Point.hh"

#include "MeshUtils.hh"

namespace Amanzi {
namespace AmanziMesh {

//
// Typedefs
//
using Entity_ID = int;
using Entity_GID = int;
using Set_ID = int;
using size_type = Kokkos::MeshView<int*, Kokkos::DefaultHostExecutionSpace>::size_type;

//
// Views are on host or device
//
template<typename T> using View_type = Kokkos::MeshView<T*, Kokkos::DefaultHostExecutionSpace>;
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


using Map_type = Epetra_Map;
using Map_ptr_type = Teuchos::RCP<Map_type>;
using Import_type = Epetra_Import;
using Import_ptr_type = Teuchos::RCP<Import_type>;


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


template<MemSpace_kind MEM>
struct MeshCache;

template<MemSpace_kind MEM = MemSpace_kind::HOST, AccessPattern_kind AP = AccessPattern_kind::DEFAULT> 
struct Getter {
  template<typename DATA, typename MF, typename FF, typename CF> 
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get(bool cached, DATA& d, MF& mf, FF&& f, CF&& c, const Entity_ID i){
      using type_t = typename DATA::t_dev::traits::value_type; 
      // To avoid the cast to non-reference 
      if (cached) return static_cast<type_t>(view<MEM>(d)(i));
      if constexpr(MEM == MemSpace_kind::HOST){
        if constexpr (!std::is_same_v<FF,decltype(nullptr)>)
          if(mf.get())
            return f(i);
      }
      if constexpr (!std::is_same_v<CF,decltype(nullptr)>)
        return c(i);
      assert(false); 
      return type_t{}; 
  }
}; // Getter

template<MemSpace_kind MEM> 
struct Getter<MEM,AccessPattern_kind::CACHE> {
  template<typename DATA, typename MF, typename FF, typename CF> 
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get(bool cached, DATA& d, MF&, FF&&, CF&&, const Entity_ID i){
      assert(cached);
      return view<MEM>(d)(i);
  }
}; // Getter

template<MemSpace_kind MEM> 
struct Getter<MEM,AccessPattern_kind::FRAMEWORK> {
  template<typename DATA, typename MF, typename FF, typename CF> 
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get(bool, DATA&, MF& mf, FF&& f, CF&&, const Entity_ID i){
      static_assert(!std::is_same<FF,decltype(nullptr)>::value); 
      static_assert(MEM == MemSpace_kind::HOST);
      assert(mf.get()); 
      return f(i);
  }
}; // Getter


template<MemSpace_kind MEM> 
struct Getter<MEM,AccessPattern_kind::COMPUTE> {
  template<typename DATA, typename MF, typename FF, typename CF> 
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get(bool, DATA&, MF&, FF&&, CF&& c, const Entity_ID i){
    static_assert(!std::is_same<CF,decltype(nullptr)>::value); 
    return c(i); 
  }
}; // Getter


// Getters for raggedViews
template<MemSpace_kind MEM = MemSpace_kind::HOST, AccessPattern_kind AP = AccessPattern_kind::DEFAULT>
struct RaggedGetter{ 
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto)
  get (bool cached, DATA& d, MF& mf, FF&& f, CFD&& cd, CF&& c, const Entity_ID n) {
    using view_t = typename decltype(d.template getRow<MEM>(n))::const_type;
    if (cached) {
      return view_t(d.template getRow<MEM>(n));  
    }
    if constexpr(MEM == MemSpace_kind::HOST){

      if constexpr (!std::is_same<FF,decltype(nullptr)>::value){
        if(mf.get()){
          view_t v = f(n); 
          return v; 
        }
      }
      if constexpr (!std::is_same<CF,decltype(nullptr)>::value){
        view_t v = c(c);
        return v; 
      }
    } else {
      if constexpr (!std::is_same<CF,decltype(nullptr)>::value) {
        view_t v = cd(c); 
        return v;
      }
    }
    assert(false && "No access to cache/framework/compute available"); 
    return view_t{}; 
  }
};

template<MemSpace_kind MEM>
struct RaggedGetter<MEM,AccessPattern_kind::CACHE>{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get (bool cached, DATA& d, MF&, FF&&, CFD&&, CF&&, const Entity_ID n) { 
    assert(cached);
    return d.template getRow<MEM>(n); 
  }
};

template<MemSpace_kind MEM>
struct RaggedGetter<MEM,AccessPattern_kind::FRAMEWORK>{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get (bool, DATA&, MF& mf, FF&& f, CFD&&, CF&&, const Entity_ID n) { 
    static_assert(!std::is_same<FF,decltype(nullptr)>::value); 
    static_assert(MEM == MemSpace_kind::HOST);
    assert(mf.get()); 
    return f(n);
  }
};

template<MemSpace_kind MEM>
struct RaggedGetter<MEM,AccessPattern_kind::COMPUTE>{
  template<typename DATA, typename MF, typename FF, typename CFD, typename CF>
  static KOKKOS_INLINE_FUNCTION decltype(auto) 
  get (bool, DATA&, MF&, FF&&, CFD&& cd, CF&& c, const Entity_ID n) { 
    if constexpr(MEM == MemSpace_kind::HOST){
      static_assert(!std::is_same<CF,decltype(nullptr)>::value); 
      return c(n); 
    }else{
      static_assert(MEM==MemSpace_kind::DEVICE); 
      static_assert(!std::is_same<CFD,decltype(nullptr)>::value); 
      return cd(n); 
    }
  }
};

using MeshSets = std::map<std::tuple<std::string,Entity_kind,Parallel_kind>, Entity_ID_DualView>;
using MeshSetVolumeFractions = std::map<std::tuple<std::string,Entity_kind,Parallel_kind>, Double_DualView>;


}  // namespace AmanziMesh
}  // namespace Amanzi


