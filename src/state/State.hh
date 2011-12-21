/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the State.  State is a simple data-manager, allowing PKs to
require, read, and write various fields.  Provides some data protection by
providing both const and non-const fields to PKs.  Provides some
initialization capability -- this is where all independent variables can be
initialized (as independent variables are owned by state, not by any PK).
------------------------------------------------------------------------- */

#ifndef __STATE_HH__
#define __STATE_HH__

#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Mesh.hh"
#include "Vis.hpp"
#include "Field.hh"

namespace Amanzi {

typedef enum { COMPLETE, UPDATING } status_type;

class State {

public:
  // types
  typedef std::vector< Teuchos::RCP<Field> > Fields;
  typedef std::map<std::string, Fields::size_type> FieldNameMap;

  State(Teuchos::RCP<Amanzi::AmanziMesh::Mesh>&);

  State(Teuchos::ParameterList&,
	 Teuchos::RCP<Amanzi::AmanziMesh::Mesh>&);

  State(const State&); // copy constructor
  State& operator=(const State&); // assignment

  ~State();

  // initialize values over blocks
  void initialize();
  bool check_all_initialized();

  // Add things to the state.  Location is one of the AmanziMesh::Entity_kind
  // enum, own indicates whether the PK requiring the field would like write
  // access to the field (i.e. the field is a dependent variable for the PK),
  // and num_dofs indicates the number of required vectors.
  //
  // Note that multiple PKs may require a field, but only one may own it.
  void require_field(std::string fieldname, int location,
                     std::string owner="state", int num_dofs=1);

  // -- access methods --
  // This access method should be used by PKs who don't own the field, i.e.
  // flow accessing a temperature field if an energy PK owns the temperature field.
  // This ensures a PK cannot mistakenly alter data it doesn't own.
  Teuchos::RCP<const Epetra_MultiVector> get_field(std::string fieldname) const;

  // This access method should be used by PKs who own the field.
  Teuchos::RCP<Epetra_MultiVector> get_field(std::string fieldname,
                                             std::string pk_name);

  // Access to the full field instance, not just the data.
  Teuchos::RCP<Field> get_field_record(std::string fieldname);

  // CRUFT, fix me.
  Teuchos::RCP<double> get_density() { return density_; }
  Teuchos::RCP<double> get_viscosity() { return viscosity_; }
  Teuchos::RCP<double*> get_gravity() { return gravity_; }

  Amanzi::AmanziMesh::Mesh& get_mesh() { return *mesh_maps_; }
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> get_mesh_maps() { return mesh_maps_; };

  double get_time () const { return time_; }
  int get_cycle () const { return cycle_; }

  // modify methods
  void set_field(std::string fieldname, std::string pk_name,
                 const Epetra_Vector&);
  void set_field(std::string fieldname, std::string pk_name,
                 const Epetra_MultiVector&);
  void set_field(std::string fieldname, std::string pk_name, const double* u);
  void set_field(std::string fieldname, std::string pk_name, double u);

  // modify methods only valid for cell-based fields
  void set_field(std::string fieldname, std::string pk_name,
                 const double* u, int mesh_block_id);
  void set_field(std::string fieldname, std::string pk_name,
                 double u, int mesh_block_id);

  // modify methods only valid for face-based fields
  void set_vector_field(std::string fieldname, std::string pk_name,
                        const double* u, int mesh_block_id);

  // modify subfield_names, which are used for io naming
  void set_subfield_names(std::string fieldname,
                          const std::vector<std::string>& subfield_names);

  void set_time ( double new_time ) { time_ = new_time; }
  void advance_time(double dT) { time_ += dT; }
  void set_cycle (int new_cycle ) { cycle_ = new_cycle; }

  // CRUFT, fix me
  void set_density(double wd);
  void set_viscosity(double mu);
  void set_gravity(const double *g);

  // status methods
  const status_type get_status () const { return status_; };
  void set_status ( status_type new_status ) { status_ = new_status; }

  // vis and restart functions
  void write_vis(Amanzi::Vis& vis);
  void write_vis(Amanzi::Vis& vis, Epetra_MultiVector *auxdata, std::vector<std::string>& auxnames);

private:
  Teuchos::RCP<const Field> get_field_record(std::string fieldname) const;

  void set_cell_value_in_mesh_block(double value, Epetra_MultiVector &v,
                                    int mesh_block_id);

  void initialize_from_parameter_list();

  // field container and fieldname map from name -> container location
  FieldNameMap field_name_map_;
  Fields fields_;

  // CRUFT, fix me
  // extra data which (for the moment) is assumed constant
  // over the entire domain
  Teuchos::RCP<double*> gravity_;
  Teuchos::RCP<double> density_;
  Teuchos::RCP<double> viscosity_;

  double time_;
  int cycle_;
  status_type status_;

  // mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh_maps_;

  // parameter list
  Teuchos::ParameterList parameter_list_;
};

inline
Teuchos::RCP<Field> State::get_field_record(std::string fieldname) {
  return fields_[field_name_map_[fieldname]];
};

inline
Teuchos::RCP<const Field> State::get_field_record(std::string fieldname) const {
  return fields_[field_name_map_.find(fieldname)->second];
};

} // namespace amanzi
#endif
