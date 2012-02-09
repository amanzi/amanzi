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

#ifndef STATE_STATE_HH_
#define STATE_STATE_HH_

#include <string>
#include <vector>
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
    State(Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh_maps);

    State(Teuchos::ParameterList& parameter_list,
          Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh_maps);

    // Copy constructor, copies memory not pointers.
    explicit State(const State& s);

    // Assignment operator, copies memory not pointers.  Note this
    // implementation requires the State being copied has the same structure (in
    // terms of fields, order of fields, etc) as *this.  This really means that
    // it should be a previously-copy-constructed version of the State.  One and
    // only one State should be instantiated and populated -- all other States
    // should be copy-constructed from that initial State.
    State& operator=(const State& s);

    // initialize values over blocks
    void initialize();
    bool check_all_initialized();

    // Add things to the state.
    //
    // Note that multiple PKs may require a field, but only one may own it.
    void require_field(std::string fieldname, std::string owner="state");

    // -- access methods --
    // This access method should be used by PKs who don't own the field, i.e.
    // flow accessing a temperature field if an energy PK owns the temperature field.
    // This ensures a PK cannot mistakenly alter data it doesn't own.
    Teuchos::RCP<const CompositeVector> get_data(std::string fieldname) const;

    // This access method should be used by PKs who own the field.
    Teuchos::RCP<CompositeVector> get_data(std::string fieldname, std::string pk_name);

    // Access to the full field instance, not just the data.
    Teuchos::RCP<Field> get_field(std::string fieldname);

    // CRUFT, fix me.
    Teuchos::RCP<double> density() { return density_; }
    Teuchos::RCP<const double> density() const { return density_; }

    Teuchos::RCP<double> viscosity() { return viscosity_; }
    Teuchos::RCP<const double> viscosity() const { return viscosity_; }

    Teuchos::RCP< std::vector<double> > gravity() { return gravity_; }
    Teuchos::RCP< const std::vector<double> > gravity() const { return gravity_; }

    Teuchos::RCP<AmanziMesh::Mesh> mesh() { return mesh_; }
    Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }

    double time () const { return time_; }
    int cycle () const { return cycle_; } // what is the cycle?
    status_type status () const { return status_; };

    // modify methods
    void set_data(std::string fieldname, std::string pk_name,
                  Teuchos::RCP<CompositeVector>& data);
    void set_data(std::string fieldname, std::string pk_name,
                   const CompositeVector& data);

    void set_time ( double new_time ) { time_ = new_time; }
    void advance_time(double dT) { time_ += dT; }
    void set_cycle (int new_cycle ) { cycle_ = new_cycle; }

    // CRUFT, fix me
    void set_density(double wd);
    void set_viscosity(double mu);
    void set_gravity(const double g[3]);
    void set_gravity(const std::vector<double> g);

    // status methods

    void set_status ( status_type new_status ) {
      status_ = new_status;
    }

    // vis and restart functions
    void write_vis(Amanzi::Vis& vis);
    void write_vis(Amanzi::Vis& vis, Epetra_MultiVector *auxdata, std::vector<std::string>& auxnames);

  private:
    Teuchos::RCP<const Field> get_field_record(std::string fieldname) const;

    void set_cell_value_in_mesh_block(double value, Epetra_MultiVector &v,
            int mesh_block_id);

    void initialize_from_parameter_list();

    // field container and fieldname map from name -> container location
    std::vector< Teuchos::RCP<Field> > fields_;
    std::map<std::string, std::vector< Teuchos::RCP<Field> >::size_type> field_name_map_;

    // CRUFT, fix me
    // extra data which (for the moment) is assumed constant
    // over the entire domain
    Teuchos::RCP< std::vector<double> > gravity_;
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
