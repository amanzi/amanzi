/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the State.  State is a simple data-manager, allowing PKs to
require, read, and write various fields.  Provides some data protection by
providing both const and non-const fields to PKs.  Provides some
initialization capability -- this is where all independent variables can be
initialized (as independent variables are owned by state, not by any PK).

------------------------------------------------------------------------- */

/* TODO: (etc 12/21), ATS ticket #6
1. Yank crufty density and viscosity out of here... they may be spatially
variable.

2. Consider making Field virtual and allowing an implementation with
scalars/NumVectors() length vectors to decrease memory footprint for things
like density which may NOT be spatially variable.
*/

#include <iostream>

#include "Epetra_Vector.h"

#include "errors.hh"
#include "CompositeVector.hh"
#include "State.hh"

//namespace Amanzi {

State::State(Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh):
  mesh_(mesh) {};

State::State(Teuchos::ParameterList& state_plist, Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
  mesh_(mesh), state_plist_(state_plist) {};

// copy constructor:
// Create a new State with different data but the same values.
//
// Could get a better implementation with a CopyMode, see TransportState in
// Amanzi as an example.  I'm not sure its needed at this point, however.
State::State(const State& other) :
  mesh_(other.mesh_), state_plist_(other.state_plist_) {

  field_name_map_ = s.field_name_map_;
  fields_.resize(s.fields_.size());
  for (unsigned int lcv = 0; lcv != s.fields_.size(); ++lcv) {
    fields_[lcv] = other.fields_[lcv]->Clone();
  }

  time_ = other.time_;
  cycle_ = other.cycle_;
};

// operator=:
//  Assign a state's data from another state.  Note this
// implementation requires the State being copied has the same structure (in
// terms of fields, order of fields, etc) as *this.  This really means that
// it should be a previously-copy-constructed version of the State.  One and
// only one State should be instantiated and populated -- all other States
// should be copy-constructed from that initial State.
State& State::operator=(const State& other) {
  if (this != &other) {
    if (fields_.size() != other.fields_.size()) {
      Errors::Message message("Attempted assignment of non-compatible states.");
      Exceptions::amanzi_throw(message);
    }

    for (unsigned int lcv = 0; lcv != fields_.size(); ++lcv) {
      fields_[lcv]->set_io_restart(other.fields_[lcv]->io_restart());
      fields_[lcv]->set_io_vis(other.fields_[lcv]->io_vis());
      fields_[lcv]->set_initialized(other.fields_[lcv]->initialized());
      if (fields_[lcv]->type() == VECTOR_FIELD) {
        fields_[lcv]->set_data(fields_[lcv]->owner(), *other.fields_[lcv]->vector_data());
      } else if (fields_[lcv]->type() == CONSTANT_VECTOR) {
        fields_[lcv]->set_data(fields_[lcv]->owner(), *other.fields_[lcv]->constant_vector_data());
      } else if (fields_[lcv]->type() == CONSTANT_SCALAR) {
        fields_[lcv]->set_data(fields_[lcv]->owner(), *other.fields_[lcv]->scalar_data());
      }
    }

    time_ = other.time_;
    cycle_ = other.cycle_;
    status_ = other.status_;
  }
  return *this;
};

// Initialize data, allowing values to be specified here or in the owning PK.
// All independent variables must be initialized here.
void State::Initialize() {
  InitializeFromParameterList_();
};

// Initialize fields from the parameter list of "Constant {Fieldname}",
// including all independent variables.
void State::InitializeFromParameterList_() {
  // gravity
  
  double u[3];
  u[0] = parameter_list_.get<double>("Gravity x");
  u[1] = parameter_list_.get<double>("Gravity y");
  u[2] = parameter_list_.get<double>("Gravity z");
  set_gravity(u);

  if (parameter_list_.isParameter("Constant water density")) {
    double u = parameter_list_.get<double>("Constant water density");
    std::cout << "initializing in state: density =" << u << std::endl;
    set_density(u);
  }
  if (parameter_list_.isParameter("Constant viscosity")) {
    double u = parameter_list_.get<double>("Constant viscosity");
    std::cout << "initializing in state: viscosity =" << u << std::endl;
    set_viscosity(u);
  }

  // Initialize -- State has "sudo" privleges and can write them all
  for (std::vector< Teuchos::RCP<Field> >::iterator field = fields_.begin();
       field != fields_.end(); ++field) {
    if ((*field)->get_location() == FIELD_LOCATION_CELL) {
      const std::vector<std::string> subfield_names = (*field)->get_subfield_names();
      if (subfield_names.size() == (*field)->get_num_dofs()) {
        double u[subfield_names.size()];
        bool got_them_all = true;
        for (int lcv=0; lcv != subfield_names.size(); ++lcv) {
          // attempt to pick out constant values for the field
          if (parameter_list_.isParameter("Constant "+subfield_names[lcv])) {
            u[lcv] = parameter_list_.get<double>("Constant "+subfield_names[lcv]);
            std::cout << "  got value:" << (*field)->get_fieldname() << "=" << u[lcv] << std::endl;
          } else {
            got_them_all = false;
            break;
          }
        }

        if (got_them_all) {
          std::cout << "got them all, assigning" << std::endl;
          (*field)->set_data((*field)->get_owner(), u);
          (*field)->set_initialized();
        }
      }
    }
  }

  int num_blocks = parameter_list_.get<int>("Number of mesh blocks");
  for (int nb=1; nb<=num_blocks; nb++) {
    std::stringstream pname;
    pname << "Mesh block " << nb;

    Teuchos::ParameterList sublist = parameter_list_.sublist(pname.str());

    int mesh_block_id = sublist.get<int>("Mesh block ID");

    for (std::vector< Teuchos::RCP<Field> >::iterator field = fields_.begin();
         field != fields_.end(); ++field) {
      std::cout << "checking for ICs for " << (*field)->get_fieldname() << std::endl;
      if ((*field)->get_location() == FIELD_LOCATION_CELL) {
        std::cout << "  checking for cell ICs for " << (*field)->get_fieldname() << std::endl;
        const std::vector<std::string> subfield_names = (*field)->get_subfield_names();
        if (subfield_names.size() == (*field)->get_num_dofs()) {
          std::cout << "  got names for " << (*field)->get_fieldname() << std::endl;
          double u[subfield_names.size()];
          bool got_them_all = true;
          for (int lcv=0; lcv != subfield_names.size(); ++lcv) {
            std::cout << subfield_names[lcv] << std::endl;
            // attempt to pick out constant values for the field
            if (sublist.isParameter("Constant "+subfield_names[lcv])) {
              u[lcv] = sublist.get<double>("Constant "+subfield_names[lcv]);
              std::cout << "  got value " << u[lcv] << std::endl;
            } else {
              got_them_all = false;
              break;
            }
          }

          if (got_them_all) {
            (*field)->set_data((*field)->get_owner(), u, mesh_block_id);
            (*field)->set_initialized();
            std::cout << "got them all, assigning" << std::endl;
          }
        }
      } else if ((*field)->get_location() == FIELD_LOCATION_FACE) {
        std::cout << "  checking for face ICs for " << (*field)->get_fieldname() << std::endl;
        if (1 == (*field)->get_num_dofs()) {
          // a vector field, try to get the 3 components
          std::string fieldname = (*field)->get_fieldname();
          if (sublist.isParameter("Constant "+fieldname+" x")) {
            double u[3];
            u[0] = sublist.get<double>("Constant "+fieldname+" x");
            u[1] = sublist.get<double>("Constant "+fieldname+" y");
            u[2] = sublist.get<double>("Constant "+fieldname+" z");

            std::cout << "initializing in state:" << fieldname << "=" << u << std::endl;
            (*field)->set_vector_data((*field)->get_owner(), u, mesh_block_id);
            (*field)->set_initialized();
          }
        }
      }
    } // loop over fields
  } // loop over blocks
};

// Make sure all fields have gotten their IC, either from State or the owning PK.
bool State::CheckAllInitialized() {
  for (std::vector< Teuchos::RCP<Field> >::iterator field = fields_.begin();
       field != fields_.end(); ++field) {
    if (!(*field)->initialized()) return false;
  }
  return true;
};

// Method for checking if the field exists or is ok to create
Teuchos::RCP<Field> State::CheckMayCreateOrOwn_or_die_(std::string fieldname, FieldType type) {
  if (field_name_map_.find(fieldname) == field_name_map_.end()) {
    // field does not exist
    return Teuchos::null;
  } else if (GetRecord_(fieldname)->owner() == "state") {
    if (type == GetRecord_(fieldname)->type()) {
      // field exists, types match, and is not owned.
      return GetRecord_(fieldname);
    } else {
      // field exists and types don't match
      std::stringstream messagestream;
      messagestream << "Requested field " << fieldname << " of type "
                    << type << " already exists as type "
                    << GetRecord_(fieldname)->type();
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
  } else {
    // Field exists and is owned
    std::stringstream messagestream;
    messagestream << "Requested field " << fieldname << " already exists and is owned by "
                  << GetRecord_(fieldname)->owner();
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

// Call the generic constructor based upon type.
void State::PushBackNewField_(std::string fieldname, FieldType type, std::string owner) {
  field_name_map_[fieldname] = fields_.size();
  if (type == CONSTANT_SCALAR) {
    fields_.push_back(Teuchos::rcp(new Field_Scalar(fieldname, owner)));
  } else if (type == CONSTANT_VECTOR) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner)));
  } else if (type == VECTOR_FIELD) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner)));
  }
};

// Field factory methods
void State::Require(std::string fieldname, FieldType type, std::string owner) {
  if (owner == "state") {
    // PK does not wish to own.
    if (field_name_map_.find(fieldname) == field_name_map_.end()) {
      // Field does not yet exist; create a new one.
      PushBackNewField_(fieldname, type, owner);

    } else if (GetRecord_(fieldname)->type != type) {
      std::stringstream messagestream;
      messagestream << "Requested field " << fieldname << " on locations "
                    << location << " already exists on location "
                    << GetRecord_(fieldname)->get_location();
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    } // else the field exists, has matching type, and PK doesn't want to own
      // anyway, so all is good.
  } else {
    Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
    if (field == Teuchos::null) {
      // Field does not yet exist; create a new one.
      PushBackNewField_(fieldname, type, owner);
    } else {
      // Field exists, is owned by state, and matches in type.  Take over ownership.
      field->set_owner(owner);
    }
  }
};

// All of the remaining constructors want to own, as they are pushing data
// information as well.

// Require scalars constant over the entire mesh.
void State::RequireScalar(std::string fieldname, std::string owner,
                          Teuchos::RCP<double>& data_ptr) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_Scalar(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};

void State::RequireScalar(std::string fieldname, std::string owner, double data) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
  Teuchos::RCP<double> data_ptr = Teuchos::rcp(new double(data));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_Scalar(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};

void State::RequireConstantVector(std::string fieldname, std::string owner,
        Teuchos::RCP<Epetra_Vector>& data_ptr) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};

void State::RequireConstantVector(std::string fieldname, std::string owner,
        const Epetra_Vector& data) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
  Teuchos::RCP<Epetra_Vector> data_ptr = Teuchos::rcp(new Epetra_Vector(data));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};

void State::RequireConstantVector(std::string fieldname, std::string owner, int dimension) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);

  // create the vector
  Epetra_LocalMap map(dimension, 0, mesh_->Comm());
  Teuchos::RCP<Epetra_Vector> data_ptr = Teuchos::rcp(new Epetra_Vector(map, false));

  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         Teuchos::RCP<CompositeVector>& data_ptr) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         const CompositeVector& data) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, type);
  Teuchos::RCP<CompositeVector> data_ptr = Teuchos::rcp(new CompositeVector(data));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->set_data(data_ptr);
  }
};


////////// DONE THROUGH HERE ///////////

  void State::RequireField(std::string fieldname, std::string owner, std::string name,
                     AmanziMesh::Entity_kind location int num_dofs=1, bool ghosted=true);
  void State::RequireField(std::string fieldname, std::string owner,
                     std::vector<std::string> names,
                     std::vector<AmanziMesh::Entity_kind> locations,
                     int num_dofs=1, bool ghosted=true);
  void RequireField(std::string fieldname, std::string owner,
                     std::vector<std::string> names,
                     std::vector<AmanziMesh::Entity_kind> locations,
                     std::vector<int> num_dofs, bool ghosted=true);





Teuchos::RCP<const Epetra_MultiVector> State::get_field(std::string fieldname) const {
  return GetRecord_(fieldname)->get_data();
};

Teuchos::RCP<Epetra_MultiVector> State::get_field(std::string fieldname,
                                                  std::string pk_name) {
  return GetRecord_(fieldname)->get_data(pk_name);
};

void State::set_field_pointer(std::string fieldname, std::string pk_name,
        Teuchos::RCP<Epetra_MultiVector>& data) {
  GetRecord_(fieldname)->set_data_pointer(pk_name, data);
};

void State::set_field(std::string fieldname, std::string pk_name,
                      const Epetra_MultiVector& data) {
  GetRecord_(fieldname)->set_data(pk_name, data);
};

void State::set_field(std::string fieldname, std::string pk_name,
                      const Epetra_Vector& data) {
  GetRecord_(fieldname)->set_data(pk_name, data);
};

void State::set_field(std::string fieldname, std::string pk_name,
                      const double* u) {
  GetRecord_(fieldname)->set_data(pk_name, u);
};

void State::set_field(std::string fieldname, std::string pk_name,
                      double u) {
  GetRecord_(fieldname)->set_data(pk_name, u);
};

void State::set_field(std::string fieldname, std::string pk_name,
                      const double* u, int mesh_block_id) {
  GetRecord_(fieldname)->set_data(pk_name, u, mesh_block_id);
};

void State::set_field(std::string fieldname, std::string pk_name,
                      double u, int mesh_block_id) {
  GetRecord_(fieldname)->set_data(pk_name, u, mesh_block_id);
};

void State::set_vector_field(std::string fieldname, std::string pk_name,
                             const double* u, int mesh_block_id) {
  GetRecord_(fieldname)->set_vector_data(pk_name, u, mesh_block_id);
};

void State::set_subfield_names(std::string fieldname,
                               const std::vector<std::string>& subfield_names) {
  GetRecord_(fieldname)->set_subfield_names(subfield_names);
};

void State::set_density(double wd ) {
  *density_ = wd;
};

void State::set_viscosity(double mu) {
  *viscosity_ = mu;
};

void State::set_gravity(const double g[3]) {
  (*gravity_)[0] = g[0];
  (*gravity_)[1] = g[1];
  (*gravity_)[2] = g[2];
};

void State::set_gravity(const std::vector<double> g) {
  *gravity_ = g;
};

void State::write_vis(Amanzi::Vis& vis) {
  if (vis.dump_requested(get_cycle()) && !vis.is_disabled()) {
    // create the new time step...
    vis.create_timestep(get_time(),get_cycle());

    // dump all the state vectors into the file
    for (std::vector< Teuchos::RCP<Field> >::iterator field = fields_.begin();
         field != fields_.end(); ++field) {
      if ((*field)->io_vis()) {
        const std::vector<std::string> subfield_names = (*field)->get_subfield_names();
        vis.write_vector(*(*field)->get_data(), subfield_names);
      }
    }
  }
};
} // namespace amanzi
