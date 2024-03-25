/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#pragma once

#include "Teuchos_RCP.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>

namespace Amanzi {
namespace AmanziMesh {

namespace Impl {

// base class for MeshAudit mixins
template <class Mesh_type>
class MeshAudit_Base {
 protected:
  using cEntity_ID_View = typename Mesh_type::cEntity_ID_View;
  using Entity_ID_View = typename Mesh_type::Entity_ID_View;

 public:
  MeshAudit_Base(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os);

 protected:
  // helpers for doing the testing
  KOKKOS_INLINE_FUNCTION
  bool areDistinctValues_(const cEntity_ID_View& list) const;

  template <class View1_type, class View2_type>
  KOKKOS_INLINE_FUNCTION
  int isSameFace_(const View1_type& l1,
                  const View2_type& l2) const;

  bool globalAny_(bool) const;

protected:
  Teuchos::RCP<const Mesh_type> mesh_;
  Comm_ptr_type comm_;

  const int nnodes_all_;
  const int nfaces_all_;
  const int ncells_all_;
  const int nnodes_owned_;
  const int nfaces_owned_;
  const int ncells_owned_;

  std::ostream& os_;
};


template <class Mesh_type>
class MeshAudit_Geometry : public MeshAudit_Base<Mesh_type> {
  using cEntity_ID_View = typename Mesh_type::cEntity_ID_View;
  using Entity_ID_View = typename Mesh_type::Entity_ID_View;

 public:
  MeshAudit_Geometry(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os = std::cout)
    : MeshAudit_Base<Mesh_type>(mesh, os)
  {}

  // The individual tests are also available.  While the tests are all formally
  // independent, there is an implicit order dependence of the tests in that a
  // test may assume certain mesh data has been verified, and that verification
  // is done by other tests.
  bool checkEntityCounts() const;
  bool checkCellToNodes() const;
  bool checkCellToFaces() const;
  bool checkFaceToNodes() const;
  bool checkCellToFaceDirs() const;
  bool checkNodeRefsByCells() const;
  bool checkFaceRefsByCells() const;
  bool checkNodeRefsByFaces() const;
  bool checkCellDegeneracy() const;
  bool checkCellGeometry() const;
  bool checkCellToFacesToNodes() const;
  bool checkNodeToCoordinates() const;
  bool checkCellToCoordinates() const;
  bool checkFaceToCoordinates() const;
  bool checkFaceCellAdjacencyConsistency() const;
  bool checkFaceNormalReltoCell() const;
  bool checkFaceNormalOrientation() const;

 protected:
  using MeshAudit_Base<Mesh_type>::mesh_;
  using MeshAudit_Base<Mesh_type>::os_;
  using MeshAudit_Base<Mesh_type>::globalAny_;
  using MeshAudit_Base<Mesh_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type>::nnodes_owned_;
};


template <class Mesh_type>
class MeshAudit_Maps : public MeshAudit_Geometry<Mesh_type> {
 protected:
  using cEntity_ID_View = typename Mesh_type::cEntity_ID_View;
  using Entity_ID_View = typename Mesh_type::Entity_ID_View;

 public:
  MeshAudit_Maps(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os = std::cout)
    : MeshAudit_Geometry<Mesh_type>(mesh, os)
  {}

  // This is the main method.
  void create_test_dependencies();

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


  using MeshAudit_Base<Mesh_type>::mesh_;
  using MeshAudit_Base<Mesh_type>::os_;
  using MeshAudit_Base<Mesh_type>::globalAny_;
  using MeshAudit_Base<Mesh_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type>::nnodes_owned_;
};


template <class Mesh_type>
class MeshAudit_Sets : public MeshAudit_Maps<Mesh_type> {
 protected:
  using cEntity_ID_View = typename Mesh_type::cEntity_ID_View;
  using Entity_ID_View = typename Mesh_type::Entity_ID_View;

 public:
  MeshAudit_Sets(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os = std::cout)
    : MeshAudit_Maps<Mesh_type>(mesh, os)
  {}

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
  using MeshAudit_Base<Mesh_type>::mesh_;
  using MeshAudit_Base<Mesh_type>::os_;
  using MeshAudit_Base<Mesh_type>::globalAny_;
  using MeshAudit_Base<Mesh_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type>::nnodes_all_;
  using MeshAudit_Base<Mesh_type>::ncells_owned_;
  using MeshAudit_Base<Mesh_type>::nfaces_owned_;
  using MeshAudit_Base<Mesh_type>::nnodes_owned_;
};


template <typename MeshAudit_type>
struct AuditGraph {
  // the tests are stored in a graph, ensuring that if a test fails, tests that
  // depend upon the capability in that test are not run.

  // method pointer
  typedef bool (MeshAudit_type::*Test)() const;

  // The vertex type for the test dependency graph.
  struct Vertex {
    Vertex() : run(true) {}
    std::string name;
    mutable bool run;
    Test test;
  };

  // The graph
  using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, Vertex>;
  Graph g;

  // things to do on the graph
  struct mark_do_not_run : public boost::bfs_visitor<> {
    template <class Vertex, class Graph>
    void discover_vertex(Vertex v, Graph& gr)
    {
      gr[v].run = false;
    }
  };

  // find a vertex by name
  typename Graph::vertex_iterator findVertex(const std::string& name)
  {
    typename Graph::vertex_iterator vi, vi_end;
    for (std::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi) {
      if (g[*vi].name == name) { return vi; }
    }
    return vi_end;
  }

  // add vertex
  typename Graph::vertex_descriptor addVertex() { return boost::add_vertex(g); }

  Vertex& getVertex(const typename Graph::vertex_descriptor& desc) { return g[desc]; }

  // add edge
  void addEdge(const std::string& vert_name_out, typename Graph::vertex_descriptor& vert_in)
  {
    auto vert_out = findVertex(vert_name_out);
    boost::add_edge(*vert_out, vert_in, g);
  }

  void
  addEdge(typename Graph::vertex_descriptor& vert_out, typename Graph::vertex_descriptor& vert_in)
  {
    boost::add_edge(vert_out, vert_in, g);
  }

  // Verify runs all the tests in an order that respects the dependencies
  // between tests.  If a particular test fails, all other tests that have it
  // as a pre-requisite are skipped.  It is important that each test return
  // a collective fail/pass result in parallel, so that all processes proceed
  // through the tests in lockstep.
  int verify(MeshAudit_type& audit, std::ostream& os) const
  {
    int status = 0;
    typedef typename Graph::vertex_descriptor GraphVertex;
    std::list<GraphVertex> run_order;
    boost::topological_sort(g, std::front_inserter(run_order));

    mark_do_not_run vis;
    for (auto itr = run_order.begin(); itr != run_order.end(); ++itr) {
      if (g[*itr].run) {
        os << "Checking " << g[*itr].name << " ..." << std::endl;
        if ((audit.*(g[*itr].test))()) {
          status = 1;
          os << "  Test failed!" << std::endl;
          boost::breadth_first_search(g, *itr, visitor(vis));
        }
      } else {
        os << "Skipping " << g[*itr].name << " check because of previous failures." << std::endl;
      }
    }
    return status;
  }
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
createTestDependencies_Base(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Geometry(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Maps(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Sets(MeshAudit_type& audit, AuditGraph<MeshAudit_type>& graph);

// Call these instead
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Base<Mesh_type>& audit,
                       AuditGraph<MeshAudit_Base<Mesh_type>>& graph)
{
  createTestDependencies_Base(audit, graph);
}
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Geometry<Mesh_type>& audit,
                       AuditGraph<MeshAudit_Geometry<Mesh_type>>& graph)
{
  createTestDependencies_Geometry(audit, graph);
}
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Maps<Mesh_type>& audit,
                       AuditGraph<MeshAudit_Maps<Mesh_type>>& graph)
{
  createTestDependencies_Maps(audit, graph);
}
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Sets<Mesh_type>& audit,
                       AuditGraph<MeshAudit_Sets<Mesh_type>>& graph)
{
  createTestDependencies_Sets(audit, graph);
}

// helper functions
template <class View_type>
bool checkErrorList(const View_type& view, const std::string& msg, std::ostream& os);

template <class View_type>
void printErrorList(const View_type& view, const std::string& msg, int count, std::ostream& os);

inline bool globalAny(const Comm_type& comm, bool value) {
  int lval = value, gval;
  Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &lval, &gval);
  return (bool) gval;
}

} // namespace Impl


template <class Mesh_type, template <typename T> class MeshAudit_type>
class MeshAudit_ {
 public:
  MeshAudit_(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os = std::cout)
    : os_(os), audit_(mesh, os), graph_()
  {
    Impl::createTestDependencies(audit_, graph_);
  }

  int Verify() { return graph_.verify(audit_, os_); }

 private:
  std::ostream& os_;
  MeshAudit_type<Mesh_type> audit_;
  Impl::AuditGraph<MeshAudit_type<Mesh_type>> graph_;
};


} // namespace AmanziMesh
} // namespace Amanzi
