#ifndef __MESHAUDIT__
#define __MESHAUDIT__

#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"

#include "Mesh.hh"

class MeshAudit {
public:

  MeshAudit(Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh_, std::ostream& os=std::cout);

  // This is the main method.
  int Verify() const;
  
  // The individual tests are also available.  While the tests are all formally
  // independent, there is an implicit order dependence of the tests in that a
  // test may assume certain mesh data has been verified, and that verification
  // is done by other tests.

  int check_entity_counts() const;
  int check_cell_to_nodes_refs() const;
  int check_cell_to_faces_refs() const;
  int check_face_to_nodes_refs() const;
  int check_cell_to_nodes_consistency() const;
  int check_cell_to_faces_consistency() const;
  int check_face_to_nodes_consistency() const;
  int check_cell_to_face_dirs_basic() const;
  int check_cell_degeneracy() const;
  int check_cell_to_faces() const;
  int check_node_to_coordinates() const;
  int check_cell_to_coordinates() const;
  int check_face_to_coordinates() const;
  int check_cell_geometry() const;
  int check_node_maps() const;
  int check_face_maps() const;
  int check_cell_maps() const;
  int check_node_to_coordinates_ghost_data() const;
  int check_face_to_nodes_ghost_data() const;
  int check_cell_to_nodes_ghost_data() const;
  int check_cell_to_faces_ghost_data() const;
  
  int check_node_sets() const;
  int check_face_sets() const;
  int check_cell_sets() const;
  
  int check_node_partition() const;
  int check_face_partition() const;

private:

  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  const Epetra_Comm &comm;
  const int MyPID;
  const int nnode;
  const int nface;
  const int ncell;

  std::ostream& os;
  unsigned int MAX_OUT;

  bool distinct_values(const Amanzi::AmanziMesh::Entity_ID_List& list) const;
  void write_list(const std::vector<unsigned int>&, unsigned int) const;
  int same_face(const Amanzi::AmanziMesh::Entity_ID_List, const Amanzi::AmanziMesh::Entity_ID_List) const;
  int check_maps(const Epetra_Map&, const Epetra_Map&) const;
  int check_sets(Amanzi::AmanziMesh::Entity_kind, const Epetra_Map&, const Epetra_Map&) const;
  int check_set_ids(Amanzi::AmanziMesh::Entity_kind) const;
  int check_set_ids_same(Amanzi::AmanziMesh::Entity_kind) const;
  int check_valid_set_id(Amanzi::AmanziMesh::Entity_kind) const;
  int check_get_set(unsigned int, Amanzi::AmanziMesh::Entity_kind, Amanzi::AmanziMesh::Parallel_type, const Epetra_Map&) const;
  int check_used_set(unsigned int, Amanzi::AmanziMesh::Entity_kind, const Epetra_Map&, const Epetra_Map&) const;
};

#endif
