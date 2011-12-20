/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#ifndef __FIELD_HH__
#define __FIELD_HH__

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Mesh.hh"

namespace Amanzi {
class Field {

public:

  Field() {};

  Field(std::string fieldname, int location,
        Teuchos::RCP<Amanzi::AmanziMesh::Mesh> &mesh_maps,
        std::string owner="state",
        int num_dofs=1);

  // access
  std::string get_fieldname() const { return fieldname_; }
  std::string get_owner() const { return owner_; }
  int get_location() const { return location_; }
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

  void set_data(std::string pk_name, const Epetra_Vector&);
  void set_data(std::string pk_name, const Epetra_MultiVector&);
  void set_data(std::string pk_name, const double* u);
  void set_data(std::string pk_name, double u);
  void set_data(std::string pk_name, const double* u, int mesh_block_id);
  void set_data(std::string pk_name, double u, int mesh_block_id);
  void set_vector_data(std::string pk_name, const double* u);
  void set_vector_data(std::string pk_name, const double* u, int mesh_block_id);

private:
  void assert_owner_or_die(std::string pk_name) const;
  void assert_location_or_die(int location) const;
  void assert_valid_block_id_or_die(int mesh_block_id) const;

  Teuchos::RCP<Epetra_MultiVector> data_;
  std::string fieldname_;
  std::string owner_;
  std::vector<std::string> subfieldnames_;
  int location_;
  int num_dofs_;
  bool io_restart_;
  bool io_vis_;
  bool initialized_;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_;

}; // class Field
} // namespace Amanzi
#endif
