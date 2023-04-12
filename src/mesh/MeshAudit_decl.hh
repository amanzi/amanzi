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
 public:
  MeshAudit_Base(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os);

 protected:
  // helpers for doing the testing
  bool
  areDistinctValues_(const AmanziMesh::View_type<const Entity_ID, MemSpace_kind::HOST>& list) const;
  void writeList_(const Entity_ID_List&) const;
  bool globalAny_(bool) const;
  int isSameFace_(const AmanziMesh::View_type<const Entity_ID, MemSpace_kind::HOST>,
                  const AmanziMesh::View_type<const Entity_ID, MemSpace_kind::HOST>) const;

 protected:
  Teuchos::RCP<const Mesh_type> mesh_;
  Comm_ptr_type comm_;

  const int nnodes_all_;
  const int nfaces_all_;
  const int ncells_all_;

  std::ostream& os_;
  unsigned int MAX_OUT;
};


template <class Mesh_type>
class MeshAudit_Geometry : public MeshAudit_Base<Mesh_type> {
 public:
  MeshAudit_Geometry(const Teuchos::RCP<const Mesh_type>& mesh, std::ostream& os = std::cout)
    : MeshAudit_Base<Mesh_type>(mesh, os)
  {}

  // The individual tests are also available.  While the tests are all formally
  // independent, there is an implicit order dependence of the tests in that a
  // test may assume certain mesh data has been verified, and that verification
  // is done by other tests.
  bool check_entity_counts() const;
  bool check_cell_to_nodes() const;
  bool check_cell_to_faces() const;
  bool check_face_to_nodes() const;
  bool check_cell_to_face_dirs() const;
  bool check_node_refs_by_cells() const;
  bool check_face_refs_by_cells() const;
  bool check_node_refs_by_faces() const;
  bool check_cell_degeneracy() const;
  bool check_cell_geometry() const;
  bool check_cell_to_faces_to_nodes() const;
  bool check_node_to_coordinates() const;
  bool check_cell_to_coordinates() const;
  bool check_face_to_coordinates() const;
  bool check_face_cell_adjacency_consistency() const;
  bool check_face_normal_relto_cell() const;
  bool check_face_normal_orientation() const;

 protected:
  using MeshAudit_Base<Mesh_type>::mesh_;
  using MeshAudit_Base<Mesh_type>::os_;
  using MeshAudit_Base<Mesh_type>::writeList_;
  using MeshAudit_Base<Mesh_type>::globalAny_;
  using MeshAudit_Base<Mesh_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type>::nnodes_all_;
};


template <class Mesh_type>
class MeshAudit_Maps : public MeshAudit_Geometry<Mesh_type> {
 public:
  MeshAudit_Maps(const Teuchos::RCP<const Mesh_type>& mesh_, std::ostream& os = std::cout)
    : MeshAudit_Geometry<Mesh_type>(mesh_, os)
  {}

  // This is the main method.
  void create_test_dependencies();

  bool check_node_maps() const;
  bool check_face_maps() const;
  bool check_cell_maps() const;
  bool check_node_to_coordinates_ghost_data() const;
  bool check_face_to_nodes_ghost_data() const;
  bool check_cell_to_nodes_ghost_data() const;
  bool check_cell_to_faces_ghost_data() const;

  bool check_maps(const Epetra_Map&, const Epetra_Map&) const;
  bool check_node_partition() const;
  bool check_face_partition() const;

 protected:
  using MeshAudit_Base<Mesh_type>::mesh_;
  using MeshAudit_Base<Mesh_type>::os_;
  using MeshAudit_Base<Mesh_type>::writeList_;
  using MeshAudit_Base<Mesh_type>::globalAny_;
  using MeshAudit_Base<Mesh_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type>::nnodes_all_;
};


template <class Mesh_type>
class MeshAudit_Sets : public MeshAudit_Maps<Mesh_type> {
 public:
  MeshAudit_Sets(const Teuchos::RCP<const Mesh_type>& mesh_, std::ostream& os = std::cout)
    : MeshAudit_Maps<Mesh_type>(mesh_, os)
  {}

  bool check_node_set_ids() const;
  bool check_face_set_ids() const;
  bool check_cell_set_ids() const;

  bool check_valid_node_set_id() const;
  bool check_valid_face_set_id() const;
  bool check_valid_cell_set_id() const;

  bool check_node_sets() const;
  bool check_face_sets() const;
  bool check_cell_sets() const;

  bool check_get_set_ids(AmanziMesh::Entity_kind) const;
  bool check_valid_set_id(AmanziMesh::Entity_kind) const;
  bool check_sets(AmanziMesh::Entity_kind, const Epetra_Map&, const Epetra_Map&) const;
  bool check_get_set(AmanziMesh::Set_ID,
                     AmanziMesh::Entity_kind,
                     AmanziMesh::Parallel_kind,
                     const Epetra_Map&) const;
  bool check_used_set(AmanziMesh::Set_ID,
                      AmanziMesh::Entity_kind,
                      const Epetra_Map&,
                      const Epetra_Map&) const;

 protected:
  using MeshAudit_Base<Mesh_type>::mesh_;
  using MeshAudit_Base<Mesh_type>::os_;
  using MeshAudit_Base<Mesh_type>::writeList_;
  using MeshAudit_Base<Mesh_type>::globalAny_;
  using MeshAudit_Base<Mesh_type>::ncells_all_;
  using MeshAudit_Base<Mesh_type>::nfaces_all_;
  using MeshAudit_Base<Mesh_type>::nnodes_all_;
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
