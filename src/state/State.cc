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

#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "CompositeVector.hh"
#include "State.hh"

namespace Amanzi {

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

  field_name_map_ = other.field_name_map_;
  fields_.resize(other.fields_.size());
  for (unsigned int lcv = 0; lcv != other.fields_.size(); ++lcv) {
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
        fields_[lcv]->SetData(fields_[lcv]->owner(), *other.fields_[lcv]->GetFieldData());
      } else if (fields_[lcv]->type() == CONSTANT_VECTOR) {
        fields_[lcv]->SetData(fields_[lcv]->owner(), *other.fields_[lcv]->GetConstantVectorData());
      } else if (fields_[lcv]->type() == CONSTANT_SCALAR) {
        fields_[lcv]->SetData(fields_[lcv]->owner(), *other.fields_[lcv]->GetScalarData());
      }
    }

    time_ = other.time_;
    cycle_ = other.cycle_;
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
  for (std::vector< Teuchos::RCP<Field> >::iterator field = fields_.begin();
       field != fields_.end(); ++field) {
    (*field)->Initialize(state_plist_);
  }
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

    } else if (GetRecord_(fieldname)->type() != type) {
      std::stringstream messagestream;
      messagestream << "Requested field " << fieldname << " of type "
                    << type << " already exists and is of type "
                    << GetRecord_(fieldname)->type();
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
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, CONSTANT_SCALAR);
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_Scalar(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireScalar(std::string fieldname, std::string owner, double data) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, CONSTANT_SCALAR);
  Teuchos::RCP<double> data_ptr = Teuchos::rcp(new double(data));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_Scalar(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireConstantVector(std::string fieldname, std::string owner,
        Teuchos::RCP<Epetra_Vector>& data_ptr) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, CONSTANT_VECTOR);
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireConstantVector(std::string fieldname, std::string owner,
        const Epetra_Vector& data) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, CONSTANT_VECTOR);
  Teuchos::RCP<Epetra_Vector> data_ptr = Teuchos::rcp(new Epetra_Vector(data));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireConstantVector(std::string fieldname, std::string owner, int dimension) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, CONSTANT_VECTOR);

  // create the vector
  Epetra_LocalMap map(dimension, 0, *mesh_->get_comm());
  Teuchos::RCP<Epetra_Vector> data_ptr = Teuchos::rcp(new Epetra_Vector(map, false));

  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_ConstantVector(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         Teuchos::RCP<CompositeVector>& data_ptr) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, VECTOR_FIELD);
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         const CompositeVector& data) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, VECTOR_FIELD);
  Teuchos::RCP<CompositeVector> data_ptr = Teuchos::rcp(new CompositeVector(data));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         AmanziMesh::Entity_kind location, int num_dofs, bool ghosted) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, VECTOR_FIELD);
  Teuchos::RCP<CompositeVector> data_ptr = Teuchos::rcp(new CompositeVector(mesh_, fieldname,
                         location, num_dofs, ghosted));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         std::vector<std::string> names,
                         std::vector<AmanziMesh::Entity_kind> locations, int num_dofs,
                         bool ghosted) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, VECTOR_FIELD);
  Teuchos::RCP<CompositeVector> data_ptr = Teuchos::rcp(new CompositeVector(mesh_, names,
                         locations, num_dofs, ghosted));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

void State::RequireField(std::string fieldname, std::string owner,
                         std::vector<std::string> names,
                         std::vector<AmanziMesh::Entity_kind> locations,
                         std::vector<int> num_dofs, bool ghosted) {
  Teuchos::RCP<Field> field = CheckMayCreateOrOwn_or_die_(fieldname, VECTOR_FIELD);
  Teuchos::RCP<CompositeVector> data_ptr = Teuchos::rcp(new CompositeVector(mesh_, names,
                         locations, num_dofs, ghosted));
  if (field == Teuchos::null) {
    fields_.push_back(Teuchos::rcp(new Field_CV(fieldname, owner, data_ptr)));
  } else {
    field->set_owner(owner);
    field->SetData(owner, data_ptr);
  }
};

Teuchos::RCP<const double> State::GetScalarData(std::string fieldname) const {
  return GetRecord_(fieldname)->GetScalarData();
};

Teuchos::RCP<double> State::GetScalarData(std::string fieldname, std::string pk_name) {
  return GetRecord_(fieldname)->GetScalarData(pk_name);
};

Teuchos::RCP<const Epetra_Vector> State::GetConstantVectorData(std::string fieldname) const {
  return GetRecord_(fieldname)->GetConstantVectorData();
};

Teuchos::RCP<Epetra_Vector> State::GetConstantVectorData(std::string fieldname,
                         std::string pk_name) {
  return GetRecord_(fieldname)->GetConstantVectorData(pk_name);
};

Teuchos::RCP<const CompositeVector> State::GetFieldData(std::string fieldname) const {
  return GetRecord_(fieldname)->GetFieldData();
};

Teuchos::RCP<CompositeVector> State::GetFieldData(std::string fieldname, std::string pk_name) {
  return GetRecord_(fieldname)->GetFieldData(pk_name);
};

Teuchos::RCP<Field> State::GetRecord(std::string fieldname, std::string pk_name) {
  if (GetRecord_(fieldname)->owner() == pk_name) {
    return GetRecord_(fieldname);
  } else {
    std::stringstream messagestream;
    messagestream << "PK " << pk_name << " is attempting write access to field " << fieldname
                  << " which is owned by " << GetRecord_(fieldname)->owner();
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
};

Teuchos::RCP<const Field> State::GetRecord(std::string fieldname) const {
  return GetRecord_(fieldname);
};

void State::SetData(std::string fieldname, std::string pk_name, Teuchos::RCP<double>& data) {
  GetRecord_(fieldname)->SetData(pk_name, data);
};

void State::SetData(std::string fieldname, std::string pk_name, Teuchos::RCP<Epetra_Vector>& data) {
  GetRecord_(fieldname)->SetData(pk_name, data);
};

void State::SetData(std::string fieldname, std::string pk_name, Teuchos::RCP<CompositeVector>& data){
  GetRecord_(fieldname)->SetData(pk_name, data);
};

void State::SetData(std::string fieldname, std::string pk_name, const double& data) {
  GetRecord_(fieldname)->SetData(pk_name, data);
};

void State::SetData(std::string fieldname, std::string pk_name, const Epetra_Vector& data) {
  GetRecord_(fieldname)->SetData(pk_name, data);
};

void State::SetData(std::string fieldname, std::string pk_name, const CompositeVector& data) {
  GetRecord_(fieldname)->SetData(pk_name, data);
};

// void State::WriteVis(Amanzi::Vis& vis) {
//   if (vis.dump_requested(get_cycle()) && !vis.is_disabled()) {
//     // create the new time step...
//     vis.create_timestep(get_time(),get_cycle());

//     // dump all the state vectors into the file
//     for (std::vector< Teuchos::RCP<Field> >::iterator field = fields_.begin();
//          field != fields_.end(); ++field) {
//       if ((*field)->io_vis()) {
//         const std::vector<std::string> subfield_names = (*field)->get_subfield_names();
//         vis.write_vector(*(*field)->get_data(), subfield_names);
//       }
//     }
//   }
// };
} // namespace amanzi
