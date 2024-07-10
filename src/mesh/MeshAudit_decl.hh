/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

/*
Developer note: This is a bit odd.  In MeshAudit tests, we will write Kokkos
parallel operations on eihter host or device, depending upon the
MemSpace_kind of MeshOrCache_type.  In either case, we will write
KOKKOS_LAMBDA which always captures its closure variables by VALUE.  If we
are on DEVICE, then MeshOrCache_type is MeshCache, which can be captured by
value and is intended to be.  If we are on HOST, then we have either a Mesh,
which we don't really need to capture by value because by reference would be
better, or a MeshFramework, which actually CAN'T be captured by value
because it is a virtual base class, and so does not have a copy constructor
(or only a virtual Clone instead).  So in those cases, MeshOrCache_type is
actually passed a POINTER to a Mesh or MeshFramework, and the pointer is
captured by value (valid on host).

So, when creating a MeshAudit_Base or similar type, create:
  1. MeshAudit_Base<MeshFramework, MeshFramework*>
  2. MeshAudit_Base<Mesh, Mesh*>
  3. MeshAudit_Base<Mesh, MeshCache>  (note no pointer)

Then, all tests capture the value-or-pointer, and must be either used
directly or dereferenced before being used.  The above get() function does
exactly this.  So, inside all KOKKOS_LAMBDAs, this object should be placed
inside get(), e.g.:

  get(m).getCellCentroid();

*/


#pragma once

#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace AmanziMesh {

namespace Impl {

template<typename T>
KOKKOS_INLINE_FUNCTION
auto&
get(const T& t) {
  if constexpr(std::is_pointer<T>::value) return *t;
  else return t;
}

// helpers for doing the testing
template<class View_type>
KOKKOS_INLINE_FUNCTION bool areDistinctValues(const View_type& list);

template <class View1_type, class View2_type>
KOKKOS_INLINE_FUNCTION int isSameFace(const View1_type& l1, const View2_type& l2);


// base class for MeshAudit mixins
//
template <class Mesh_type, class MeshOrCache_type>
class MeshAudit_Base {
 public:
  using cEntity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::cEntity_ID_View;
  using Entity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::Entity_ID_View;

  MeshAudit_Base(const Teuchos::RCP<const Mesh_type>& mesh,
                 const MeshOrCache_type& mesh_or_cache,
                 std::ostream& os);
  MeshAudit_Base(const MeshAudit_Base& other) = delete;

 protected:
  bool globalAny_(bool) const;

 protected:
  Teuchos::RCP<const Mesh_type> mesh_;
  const MeshOrCache_type& mc_;
  Comm_ptr_type comm_;

  const int nnodes_all_;
  const int nfaces_all_;
  const int ncells_all_;
  const int nnodes_owned_;
  const int nfaces_owned_;
  const int ncells_owned_;

  std::ostream& os_;
};


template <class Mesh_type, class MeshOrCache_type>
class MeshAudit_Topology : public MeshAudit_Base<Mesh_type, MeshOrCache_type> {
 public:
  using cEntity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::cEntity_ID_View;
  using Entity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::Entity_ID_View;

  MeshAudit_Topology(const Teuchos::RCP<const Mesh_type>& mesh,
                     const MeshOrCache_type& mesh_or_cache,
                     std::ostream& os = std::cout)
    : MeshAudit_Base<Mesh_type, MeshOrCache_type>(mesh, mesh_or_cache, os)
  {}
  MeshAudit_Topology(const MeshAudit_Topology& other) = delete;

  bool checkEntityCounts() const;
  bool checkCellToNodes() const;
  bool checkCellToFaces() const;
  bool checkFaceToNodes() const;
  bool checkCellToFaceDirs() const;
  bool checkNodeRefsByCells() const;
  bool checkFaceRefsByCells() const;
  bool checkNodeRefsByFaces() const;
  bool checkCellDegeneracy() const;
  bool checkCellToFacesToNodes() const;
  bool checkNodeToCoordinates() const;
  bool checkCellToCoordinates() const;
  bool checkFaceToCoordinates() const;
  bool checkFaceCellAdjacencyConsistency() const;

 protected:
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mesh_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mc_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::os_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::globalAny_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_owned_;
};


template <class Mesh_type, class MeshOrCache_type>
class MeshAudit_Geometry : public MeshAudit_Topology<Mesh_type, MeshOrCache_type> {
 public:
  using cEntity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::cEntity_ID_View;
  using Entity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::Entity_ID_View;

  MeshAudit_Geometry(const Teuchos::RCP<const Mesh_type>& mesh,
                     const MeshOrCache_type& mesh_or_cache,
                     std::ostream& os = std::cout)
    : MeshAudit_Topology<Mesh_type, MeshOrCache_type>(mesh, mesh_or_cache, os)
  {}
  MeshAudit_Geometry(const MeshAudit_Geometry& other) = delete;

  // The individual tests are also available.  While the tests are all formally
  // independent, there is an implicit order dependence of the tests in that a
  // test may assume certain mesh data has been verified, and that verification
  // is done by other tests.
  bool checkCellGeometry() const;
  bool checkFaceNormalReltoCell() const;
  bool checkFaceNormalOrientation() const;

 protected:
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mesh_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mc_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::os_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::globalAny_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_owned_;
};


template <class Mesh_type, class MeshOrCache_type>
class MeshAudit_Maps : public MeshAudit_Geometry<Mesh_type, MeshOrCache_type> {
 public:
  using cEntity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::cEntity_ID_View;
  using Entity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::Entity_ID_View;

  MeshAudit_Maps(const Teuchos::RCP<const Mesh_type>& mesh,
                 const MeshOrCache_type& mesh_or_cache,
                 std::ostream& os = std::cout)
    : MeshAudit_Geometry<Mesh_type, MeshOrCache_type>(mesh, mesh_or_cache, os)
  {}
  MeshAudit_Maps(const MeshAudit_Maps& other) = delete;

  bool checkEntityCountsAgainstMaps() const;
  bool checkNodeMaps() const;
  bool checkFaceMaps() const;
  bool checkCellMaps() const;
  bool checkNodeToCoordinatesGhostData() const;
  bool checkFaceToNodesGhostData() const;
  bool checkCellToNodesGhostData() const;
  bool checkCellToFacesGhostData() const;

  bool checkNodePartition() const;
  bool checkFacePartition() const;

 protected:
  bool checkMaps_(const Map_type&, const Map_type&) const;


  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mesh_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mc_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::os_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::globalAny_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_owned_;
};


template <class Mesh_type, class MeshOrCache_type>
class MeshAudit_Sets : public MeshAudit_Maps<Mesh_type, MeshOrCache_type> {
 public:
  using cEntity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::cEntity_ID_View;
  using Entity_ID_View = typename std::remove_pointer<MeshOrCache_type>::type::Entity_ID_View;

  MeshAudit_Sets(const Teuchos::RCP<const Mesh_type>& mesh,
                 const MeshOrCache_type& mesh_or_cache,
                 std::ostream& os = std::cout)
    : MeshAudit_Maps<Mesh_type, MeshOrCache_type>(mesh, mesh_or_cache, os)
  {}
  MeshAudit_Sets(const MeshAudit_Sets& other) = delete;

  bool checkNodeSetIDs() const;
  bool checkFaceSetIDs() const;
  bool checkCellSetIDs() const;

  bool checkValidNodeSetID() const;
  bool checkValidFaceSetID() const;
  bool checkValidCellSetID() const;

  bool checkNodeSets() const;
  bool checkFaceSets() const;
  bool checkCellSets() const;

 protected:
  bool checkGetSetIDs_(AmanziMesh::Entity_kind) const;
  bool checkValidSetID_(AmanziMesh::Entity_kind) const;

  bool checkSets_(AmanziMesh::Entity_kind, const Map_type&, const Map_type&) const;
  bool checkGetSet_(AmanziMesh::Set_ID,
                    AmanziMesh::Entity_kind,
                    AmanziMesh::Parallel_kind,
                    const Map_type&) const;
  bool checkUsedSet_(AmanziMesh::Set_ID,
                     AmanziMesh::Entity_kind,
                     const Map_type&,
                     const Map_type&) const;

 protected:
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mesh_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::mc_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::os_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::globalAny_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type, MeshOrCache_type>::nnodes_owned_;
};


template <typename MeshAudit_type>
struct AuditDirectedGraph {
  enum Status_kind : int { NONE = 0, GOOD = 1, FAIL = 2, SKIP = 3 };

  // method pointer
  typedef bool (MeshAudit_type::*Test)() const;

  // The vertex type for the test dependency graph.
  struct Vertex {
    Vertex() : status(Status_kind::NONE) {}
    mutable int status;
    std::string name;
    Test test;
  };

  // the tests are stored in a graph, ensuring that if a test fails, tests that
  // depend upon the capability in that test are not run.

  // vertex operations
  int findVertex(const std::string& name) const
  {
    for (int i = 0; i < vertices_.size(); ++i) {
      if (vertices_[i].name == name) return i;
    }
    return -1;
  }

  int AddVertex(const std::string& name, const Test& test)
  {
    Vertex v;
    v.name = name;
    v.test = test;
    vertices_.push_back(v);
    return vertices_.size() - 1;
  }

  // edge operations
  void AddEdge(const std::string& name_out, int vert_in)
  {
    int vert_out = findVertex(name_out);
    edges_.push_back(std::make_pair(vert_out, vert_in));
  }

  void AddEdge(int vert_out, int vert_in) { edges_.push_back(std::make_pair(vert_out, vert_in)); }

  // graph operations
  int FindAnyRoot() const
  {
    int nv = vertices_.size();
    std::vector<int> flag(nv);

    for (int i = 0; i < nv; ++i) {
      if (vertices_[i].status == Status_kind::NONE) flag[i] = 1;
    }

    for (auto it = edges_.begin(); it != edges_.end(); ++it) {
      auto& v1 = vertices_[it->first];
      if (v1.status == Status_kind::NONE) flag[it->second] = 0;
    }

    for (int i = 0; i < nv; ++i)
      if (flag[i] == 1) return i;
    return -1;
  }

  void SkipBranches(int v, std::ostream& os) const
  {
    bool found(true);

    while (found) {
      found = false;
      for (auto it = edges_.begin(); it != edges_.end(); ++it) {
        auto& v1 = vertices_[it->first];
        auto& v2 = vertices_[it->second];
        if (v1.status == Status_kind::FAIL && v2.status == Status_kind::NONE) {
          v2.status = Status_kind::SKIP;
          os << "Skipping " << v2.name << " check because of previous failures." << std::endl;
          found = true;
        } else if (v1.status == Status_kind::SKIP && v2.status == Status_kind::NONE) {
          v2.status = Status_kind::SKIP;
          os << "Skipping " << v2.name << " check because of previous failures." << std::endl;
          found = true;
        }
      }
    }
  }

  // Verify runs all the tests in an order that respects the dependencies
  // between tests.  If a particular test fails, all other tests that have it
  // as a pre-requisite are skipped.  It is important that each test return
  // a collective fail/pass result in parallel, so that all processes proceed
  // through the tests in lockstep.
  int Verify(MeshAudit_type& audit, std::ostream& os) const
  {
    int ok(0), v;

    while ((v = FindAnyRoot()) != -1) {
      os << "Checking " << vertices_[v].name << " ..." << std::endl;
      if (!(audit.*(vertices_[v].test))()) {
        vertices_[v].status = Status_kind::GOOD;
      } else {
        os << "  Test failed!" << std::endl;
        vertices_[v].status = Status_kind::FAIL;
        SkipBranches(v, os);
        ok = 1;
      }
    }
    return ok;
  }

 private:
  // adjacency structure
  std::vector<Vertex> vertices_;
  std::vector<std::pair<int, int>> edges_;
};


//
// Note the difference here -- the first set take _any_ type audit and add the
// appropriate tests.  The second set take a specific type of audit and add
// those tests.
//
// Don't call the former -- the latter will call the former.
//
template <class MeshAudit_type>
void
createTestDependencies_Base(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Topology(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Geometry(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);


template <class MeshAudit_type>
void
createTestDependencies_Maps(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Sets(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

// Call these instead
template <class Mesh_type, class MeshOrCache_type>
void
createTestDependencies(MeshAudit_Base<Mesh_type, MeshOrCache_type>& audit,
                       AuditDirectedGraph<MeshAudit_Base<Mesh_type, MeshOrCache_type>>& graph)
{
  createTestDependencies_Base(audit, graph);
}
template <class Mesh_type, class MeshOrCache_type>
void
createTestDependencies(MeshAudit_Topology<Mesh_type, MeshOrCache_type>& audit,
                       AuditDirectedGraph<MeshAudit_Topology<Mesh_type, MeshOrCache_type>>& graph)
{
  createTestDependencies_Topology(audit, graph);
}
template <class Mesh_type, class MeshOrCache_type>
void
createTestDependencies(MeshAudit_Geometry<Mesh_type, MeshOrCache_type>& audit,
                       AuditDirectedGraph<MeshAudit_Geometry<Mesh_type, MeshOrCache_type>>& graph)
{
  createTestDependencies_Geometry(audit, graph);
}
template <class Mesh_type, class MeshOrCache_type>
void
createTestDependencies(MeshAudit_Maps<Mesh_type, MeshOrCache_type>& audit,
                       AuditDirectedGraph<MeshAudit_Maps<Mesh_type, MeshOrCache_type>>& graph)
{
  createTestDependencies_Maps(audit, graph);
}
template <class Mesh_type, class MeshOrCache_type>
void
createTestDependencies(MeshAudit_Sets<Mesh_type, MeshOrCache_type>& audit,
                       AuditDirectedGraph<MeshAudit_Sets<Mesh_type, MeshOrCache_type>>& graph)
{
  createTestDependencies_Sets(audit, graph);
}

// helper functions
template <class View_type>
bool
checkErrorList(const View_type& view, const std::string& msg, std::ostream& os);

template <class View_type>
void
printErrorList(const View_type& view, const std::string& msg, int count, std::ostream& os);

inline bool
globalAny(const Comm_type& comm, bool value)
{
  int lval = value, gval;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &lval, &gval);
  return (bool)gval;
}

} // namespace Impl


template <class Mesh_type, class MeshOrCache_type,
          template <typename T1, typename T2> class MeshAudit_type>
class MeshAudit_ {
 public:
  MeshAudit_(const Teuchos::RCP<const Mesh_type>& mesh,
             const MeshOrCache_type& mc,
             std::ostream& os = std::cout)
    : os_(os), audit_(mesh, mc, os), graph_()
  {
    Impl::createTestDependencies(audit_, graph_);
  }

  int Verify() { return graph_.Verify(audit_, os_); }

 private:
  std::ostream& os_;
  MeshAudit_type<Mesh_type, MeshOrCache_type> audit_;
  Impl::AuditDirectedGraph<MeshAudit_type<Mesh_type, MeshOrCache_type>> graph_;
};


} // namespace AmanziMesh
} // namespace Amanzi
