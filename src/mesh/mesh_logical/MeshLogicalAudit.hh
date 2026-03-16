/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef __MESHAUDIT__
#define __MESHAUDIT__

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

namespace Amanzi {
namespace AmanziMesh {

class MeshLogicalAudit {
 public:
  enum Status_kind : int { NONE = 0, GOOD = 1, FAIL = 2, SKIP = 3 };

  MeshLogicalAudit(const Teuchos::RCP<const AmanziMesh::MeshHost>& mesh_,
                   std::ostream& os = std::cout);

  // This is the main method.
  int Verify() const;

  // The individual tests are also available.  While the tests are all formally
  // independent, there is an implicit order dependence of the tests in that a
  // test may assume certain mesh data has been verified, and that verification
  // is done by other tests.

  bool check_entity_counts() const;
  bool check_cell_to_faces() const;
  bool check_face_refs_by_cells() const;
  bool check_cell_to_face_dirs() const;
  bool check_faces_cell_consistency() const;
  bool check_cell_degeneracy() const;
  bool check_cell_geometry() const;
  bool check_face_geometry() const;
  bool check_cell_face_geometry() const;
  bool check_face_maps() const;
  bool check_cell_maps() const;
  bool check_cell_to_faces_ghost_data() const;
  bool check_face_partition() const;
  bool check_cell_face_bisector_geometry() const;

 private:
  Teuchos::RCP<const AmanziMesh::MeshHost> mesh;

  Comm_ptr_type comm_;
  const int MyPID;
  const int nface;
  const int ncell;

  std::ostream& os;
  unsigned int MAX_OUT;

  bool distinct_values(const MeshHost::cEntity_ID_View& list) const;
  void write_list(const Entity_ID_List&, unsigned int) const;
  bool global_any(bool) const;
  int same_face(const MeshHost::Entity_ID_View, const MeshHost::Entity_ID_View) const;

  bool check_maps(const Epetra_Map&, const Epetra_Map&) const;
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

  // This is the vertex type for the test dependency graph.
  typedef bool (MeshLogicalAudit::*Test)() const;

  struct Vertex {
    Vertex()
      : status(Status_kind::NONE)
    {}
    std::string name;
    mutable int status;
    Test test;
  };

  // vertex operations
  int AddVertex(const std::string& name, const Test& test)
  {
    Vertex v;
    v.name = name;
    v.test = test;
    vertices_.push_back(v);
    return vertices_.size() - 1;
  }

  // edge operations
  void AddEdge(int vert_out, int vert_in) { edges_.emplace_back(vert_out, vert_in); }

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

  void create_test_dependencies();

 private:
  // adjacency structure
  std::vector<Vertex> vertices_;
  std::vector<std::pair<int, int>> edges_;
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
