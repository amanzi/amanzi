/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#pragma once

#include "Teuchos_RCP.hpp"

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
class AuditDirectedGraph {
 public:
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
createTestDependencies_Geometry(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Maps(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

template <class MeshAudit_type>
void
createTestDependencies_Sets(MeshAudit_type& audit, AuditDirectedGraph<MeshAudit_type>& graph);

// Call these instead
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Base<Mesh_type>& audit,
                       AuditDirectedGraph<MeshAudit_Base<Mesh_type>>& graph)
{
  createTestDependencies_Base(audit, graph);
}
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Geometry<Mesh_type>& audit,
                       AuditDirectedGraph<MeshAudit_Geometry<Mesh_type>>& graph)
{
  createTestDependencies_Geometry(audit, graph);
}
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Maps<Mesh_type>& audit,
                       AuditDirectedGraph<MeshAudit_Maps<Mesh_type>>& graph)
{
  createTestDependencies_Maps(audit, graph);
}
template <class Mesh_type>
void
createTestDependencies(MeshAudit_Sets<Mesh_type>& audit,
                       AuditDirectedGraph<MeshAudit_Sets<Mesh_type>>& graph)
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

  int Verify() { return graph_.Verify(audit_, os_); }

 private:
  std::ostream& os_;
  MeshAudit_type<Mesh_type> audit_;
  Impl::AuditDirectedGraph<MeshAudit_type<Mesh_type>> graph_;
};


} // namespace AmanziMesh
} // namespace Amanzi
