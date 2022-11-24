#ifndef __MESHAUDIT__
#define __MESHAUDIT__

#include "Teuchos_RCP.hpp"


#include "Mesh.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>

namespace Amanzi {

class MeshAudit {
 public:
  MeshAudit(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_, std::ostream& os = std::cout);

  // This is the main method.
  int Verify() const;

  // The individual tests are also available.  While the tests are all formally
  // independent, there is an implicit order dependence of the tests in that a
  // test may assume certain mesh data has been verified, and that verification
  // is done by other tests.

  bool check_entity_counts() const;
  bool check_cell_to_nodes() const;
  bool check_cell_to_faces() const;
  bool check_face_to_nodes() const;
  bool check_cell_to_face_dirs() const;
  bool check_cell_to_nodes_consistency() const;
  bool check_cell_to_faces_consistency() const;
  bool check_face_to_nodes_consistency() const;
  bool check_cell_to_face_dirs_consistency() const;
  bool check_node_refs_by_cells() const;
  bool check_face_refs_by_cells() const;
  bool check_node_refs_by_faces() const;
  bool check_cell_degeneracy() const;
  bool check_cell_geometry() const;
  bool check_cell_to_faces_to_nodes() const;
  bool check_node_to_coordinates() const;
  bool check_cell_to_coordinates() const;
  bool check_face_to_coordinates() const;
  bool check_node_maps() const;
  bool check_face_maps() const;
  bool check_cell_maps() const;
  bool check_node_to_coordinates_ghost_data() const;
  bool check_face_to_nodes_ghost_data() const;
  bool check_cell_to_nodes_ghost_data() const;
  bool check_cell_to_faces_ghost_data() const;

  bool check_node_set_ids() const;
  bool check_face_set_ids() const;
  bool check_cell_set_ids() const;

  bool check_valid_node_set_id() const;
  bool check_valid_face_set_id() const;
  bool check_valid_cell_set_id() const;

  bool check_node_sets() const;
  bool check_face_sets() const;
  bool check_cell_sets() const;

  bool check_node_partition() const;
  bool check_face_partition() const;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  Comm_ptr_type comm_;

  const int MyPID;
  const int nnode;
  const int nface;
  const int ncell;

  std::ostream& os;
  unsigned int MAX_OUT;

  bool distinct_values(const AmanziMesh::Entity_ID_List& list) const;
  void write_list(const AmanziMesh::Entity_ID_List&, unsigned int) const;
  bool global_any(bool) const;
  int same_face(const AmanziMesh::Entity_ID_List, const AmanziMesh::Entity_ID_List) const;

  bool check_maps(const Epetra_Map&, const Epetra_Map&) const;
  bool check_get_set_ids(AmanziMesh::Entity_kind) const;
  bool check_valid_set_id(AmanziMesh::Entity_kind) const;
  bool check_sets(AmanziMesh::Entity_kind, const Epetra_Map&, const Epetra_Map&) const;
  bool check_get_set(AmanziMesh::Set_ID,
                     AmanziMesh::Entity_kind,
                     AmanziMesh::Parallel_type,
                     const Epetra_Map&) const;
  bool check_used_set(AmanziMesh::Set_ID,
                      AmanziMesh::Entity_kind,
                      const Epetra_Map&,
                      const Epetra_Map&) const;

  // This is the vertex type for the test dependency graph.
  typedef bool (MeshAudit::*Test)() const;
  struct Vertex {
    Vertex() : run(true) {}
    std::string name;
    mutable bool run;
    Test test;
  };

  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS, Vertex> Graph;
  Graph g;

  struct mark_do_not_run : public boost::bfs_visitor<> {
    template <class Vertex, class Graph>
    void discover_vertex(Vertex v, Graph& gr)
    {
      gr[v].run = false;
    }
  };

  void create_test_dependencies();
};

} // namespace Amanzi

#endif
