#ifndef __MESHAUDIT__
#define __MESHAUDIT__

#include "Teuchos_RCP.hpp"
#include "Epetra_Comm.h"

#include "Mesh_maps_base.hh"

class MeshAudit {
public:

  MeshAudit(Teuchos::RCP<Mesh_maps_base> &mesh_, std::ostream& os=std::cout);

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
  int check_cell_topology() const;
  int check_node_maps() const;
  int check_face_maps() const;
  int check_cell_maps() const;
  int check_node_to_coordinates_ghost_data() const;
  int check_face_to_nodes_ghost_data() const;
  int check_cell_to_nodes_ghost_data() const;
  int check_cell_to_faces_ghost_data() const;

private:

  Teuchos::RCP<Mesh_maps_base> mesh;

  const Epetra_Comm *comm;
  const int MyPID;
  const int nnode;
  const int nface;
  const int ncell;

  std::ostream& os;
  unsigned int MAX_OUT;

  bool distinct_values(const std::vector<unsigned int> &list) const;
  void write_list(const std::vector<unsigned int>&, unsigned int) const;
  int same_face(const std::vector<unsigned int>, const std::vector<unsigned int>) const;
  int check_maps(const Epetra_Map&, const Epetra_Map&) const;
};

#endif
