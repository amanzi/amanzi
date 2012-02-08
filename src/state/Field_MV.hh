/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for a Field.  Field is not intended so much to hide implementation
of data as to restrict write access to data.  It freely passes out pointers to
its private data, but only passes out read-only const pointers unless you have
the secret password (AKA the name of the process kernel that owns the data).

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#ifndef STATE_FIELD_MV_HH_
#define STATE_FIELD_MV_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "Mesh.hh"

namespace Amanzi {

enum FieldLocation {
  FIELD_LOCATION_MESH, // constant on the mesh (scalar data)
  FIELD_LOCATION_CELL, // cell centers
  FIELD_LOCATION_FACE // cell faces
};

class Field {

public:
  // null constructor
  Field() {};

  // standard constructor
  Field(std::string fieldname, FieldLocation location,
        Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh_maps,
        std::string owner="state",
        int num_dofs=1);

  // copy constructor and assignment
  Field(const Field&);
  Field& operator=(const Field&);

  // access
  std::string get_fieldname() const { return fieldname_; }
  std::string get_owner() const { return owner_; }
  FieldLocation get_location() const { return location_; }
  int get_num_dofs() const { return num_dofs_; }
  const std::vector<std::string>& get_subfield_names() const { return subfieldnames_; }
  bool io_restart() const { return io_restart_; }
  bool io_vis() const { return io_vis_; }
  bool initialized() const { return initialized_; }
  Teuchos::RCP<const Epetra_MultiVector> get_data() const { return data_; }
  Teuchos::RCP<Epetra_MultiVector> get_data(std::string pk_name);
  Amanzi::AmanziMesh::Mesh& get_mesh() { return *mesh_maps_; }

  // modify
  void set_owner(std::string owner) { owner_ = owner; }
  void set_subfield_names(const std::vector<std::string>& subfieldnames) {
    subfieldnames_ = subfieldnames; }
  void set_io_restart(bool io_restart) { io_restart_ = io_restart; }
  void set_io_vis(bool io_vis) { io_vis_ = io_vis; }
  void set_initialized(bool initialized=true) { initialized_ = initialized; }

  void set_data_pointer(std::string pk_name, Teuchos::RCP<Epetra_MultiVector>& data);
  void set_data(std::string pk_name, const Epetra_Vector&);
  void set_data(std::string pk_name, const Epetra_MultiVector&);
  void set_data(std::string pk_name, const double* u);
  void set_data(std::string pk_name, double u);
  void set_data(std::string pk_name, const double* u, int mesh_block_id);
  void set_data(std::string pk_name, double u, int mesh_block_id);
  void set_vector_data(std::string pk_name, const double* u);
  void set_vector_data(std::string pk_name, const double* u, int mesh_block_id);

private:
  // consistency checking
  void assert_owner_or_die(std::string pk_name) const;
  void assert_location_or_die(FieldLocation location) const;
  void assert_valid_block_id_or_die(int mesh_block_id) const;

  Teuchos::RCP<Epetra_MultiVector> data_;
  std::string fieldname_;
  std::string owner_;
  std::vector<std::string> subfieldnames_;
  FieldLocation location_;
  int num_dofs_;
  bool io_restart_;
  bool io_vis_;
  bool initialized_;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_;

}; // class Field
} // namespace Amanzi

#endif
