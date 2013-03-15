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

#include <iostream>

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "composite_vector.hh"
#include "function-factory.hh"
#include "function.hh"
#include "field_evaluator_factory.hh"
#include "cell_volume_evaluator.hh"
#include "rank_evaluator.hh"

#include "state.hh"

namespace Amanzi {

State::State() {};

State::State(Teuchos::ParameterList& state_plist) :
    state_plist_(state_plist),
    time_(0.0),
    cycle_(0) {};

// copy constructor:
// Create a new State with different data but the same values.
//
// Could get a better implementation with a CopyMode, see TransportState in
// Amanzi as an example.  I'm not sure its needed at this point, however.
State::State(const State& other) :
    state_plist_(other.state_plist_),
    meshes_(other.meshes_),
    field_factories_(other.field_factories_),
    time_(other.time_),
    cycle_(other.cycle_) {

  for (FieldMap::const_iterator f_it=other.fields_.begin();
       f_it!=other.fields_.end(); ++f_it) {
    fields_[f_it->first] = f_it->second->Clone();
  }

  for (FieldEvaluatorMap::const_iterator fm_it=other.field_evaluators_.begin();
       fm_it!=other.field_evaluators_.end(); ++fm_it) {
    field_evaluators_[fm_it->first] = fm_it->second->Clone();
  }
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
    for (FieldMap::const_iterator f_it=other.fields_.begin();
         f_it!=other.fields_.end(); ++f_it) {
      Teuchos::RCP<Field> myfield = GetField_(f_it->first);
      if (myfield == Teuchos::null) {
        myfield = f_it->second->Clone();
        fields_[f_it->first] = myfield;
      }

      Teuchos::RCP<const Field> otherfield = f_it->second;
      myfield->set_io_checkpoint(otherfield->io_checkpoint());
      myfield->set_io_vis(otherfield->io_vis());
      myfield->set_initialized(otherfield->initialized());

      if (myfield->type() == COMPOSITE_VECTOR_FIELD) {
        myfield->SetData(*otherfield->GetFieldData());
      } else if (myfield->type() == CONSTANT_VECTOR) {
        myfield->SetData(*otherfield->GetConstantVectorData());
      } else if (myfield->type() == CONSTANT_SCALAR) {
        myfield->SetData(*otherfield->GetScalarData());
      }
    }

    for (FieldEvaluatorMap::const_iterator fm_it=other.field_evaluators_.begin();
         fm_it!=other.field_evaluators_.end(); ++fm_it) {
      Teuchos::RCP<FieldEvaluator> myfm = GetFieldEvaluator_(fm_it->first);
      if (myfm == Teuchos::null) {
        myfm = fm_it->second->Clone();
        field_evaluators_[fm_it->first] = myfm;
      }
      *myfm = *fm_it->second;
    }

    time_ = other.time_;
    cycle_ = other.cycle_;
    meshes_ = other.meshes_;
  }
  return *this;
};


// -----------------------------------------------------------------------------
// State handles mesh management.
// -----------------------------------------------------------------------------
void State::RegisterDomainMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  RegisterMesh("domain", mesh);
}


void State::RegisterMesh(Key key,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  meshes_.insert(std::make_pair(key, mesh));
};


void State::RemoveMesh(Key key) {
  meshes_.erase(key);
};


Teuchos::RCP<const AmanziMesh::Mesh> State::GetMesh(Key key) const {
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = GetMesh_(key);
  if (mesh == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Mesh " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return mesh;
};


Teuchos::RCP<const AmanziMesh::Mesh> State::GetMesh_(Key key) const {
  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second;
  } else {
    return Teuchos::null;
  }
};


// -----------------------------------------------------------------------------
// State handles data evaluation.
// -----------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator>
State::RequireFieldEvaluator(Key key) {
  Teuchos::RCP<FieldEvaluator> evaluator = GetFieldEvaluator_(key);

  // See if the key is provided by another existing evaluator.
  if (evaluator == Teuchos::null) {
    for (evaluator_iterator f_it = field_evaluator_begin();
         f_it != field_evaluator_end(); ++f_it) {
      if (f_it->second->ProvidesKey(key)) {
        evaluator = f_it->second;
        SetFieldEvaluator(key, evaluator);
        break;
      }
    }
  }

  // Get the evaluator from state's Plist
  if (evaluator == Teuchos::null) {
    // -- Get the Field Evaluator plist
    Teuchos::ParameterList fm_plist;
    if (state_plist_.isSublist("field evaluators")) {
      fm_plist = state_plist_.sublist("field evaluators");
    }

    if (fm_plist.isSublist(key)) {
      // -- Get this evaluator's plist.
      Teuchos::ParameterList sublist = fm_plist.sublist(key);
      sublist.set<Key>("evaluator name", key);

      // -- Get the model plist.
      Teuchos::ParameterList model_plist;
      if (state_plist_.isSublist("model parameters")) {
        model_plist = state_plist_.sublist("model parameters");
      }

      // -- Insert any model parameters.
      if (sublist.isParameter("model parameters")) {
        std::string modelname = sublist.get<std::string>("model parameters");
        Teuchos::ParameterList modellist = GetModelParameters(modelname);
        std::string modeltype = modellist.get<std::string>("model type");
        sublist.set(modeltype, modellist);
      } else if (sublist.isParameter("models parameters")) {
        Teuchos::Array<std::string> modelnames =
            sublist.get<Teuchos::Array<std::string> >("models parameters");
        for (Teuchos::Array<std::string>::const_iterator modelname=modelnames.begin();
             modelname!=modelnames.end(); ++modelname) {
          Teuchos::ParameterList modellist = GetModelParameters(*modelname);
          std::string modeltype = modellist.get<std::string>("model type");
          sublist.set(modeltype, modellist);
        }
      }

      // -- Create and set the evaluator.
      FieldEvaluatorFactory evaluator_factory;
      evaluator = evaluator_factory.createFieldEvaluator(sublist);
      SetFieldEvaluator(key, evaluator);
    }
  }

  // Try a cell_volume.
  if (evaluator == Teuchos::null) {
    Key cell_vol("cell_volume");
    if (key.length() >= cell_vol.length() &&
        (0 == key.compare(key.length()-cell_vol.length(), cell_vol.length(), cell_vol))) {
      Teuchos::ParameterList model_plist = state_plist_.sublist("model parameters");
      Teuchos::ParameterList plist = model_plist.sublist(key);
      plist.set("evaluator name", key);
      evaluator = Teuchos::rcp(new CellVolumeEvaluator(plist));
      SetFieldEvaluator(key, evaluator);
    }
  }


  if (evaluator == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Model for field " << key << " cannot be created in State.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return evaluator;
}


Teuchos::RCP<FieldEvaluator>
State::RequireFieldEvaluator(Key key, Teuchos::ParameterList& plist) {
  Teuchos::RCP<FieldEvaluator> evaluator = GetFieldEvaluator_(key);

  // See if the key is provided by another existing evaluator.
  if (evaluator == Teuchos::null) {
    for (evaluator_iterator f_it = field_evaluator_begin();
         f_it != field_evaluator_end(); ++f_it) {
      if (f_it->second->ProvidesKey(key)) {
        evaluator = f_it->second;
        SetFieldEvaluator(key, evaluator);
      }
    }
  }

  // Create a new evaluator.
  if (evaluator == Teuchos::null) {
    // -- Create and set the evaluator.
    FieldEvaluatorFactory evaluator_factory;
    evaluator = evaluator_factory.createFieldEvaluator(plist);
    SetFieldEvaluator(key, evaluator);
  }
  return evaluator;
}


Teuchos::RCP<FieldEvaluator> State::GetFieldEvaluator(Key key) {
  Teuchos::RCP<FieldEvaluator> evaluator = GetFieldEvaluator_(key);
  if (evaluator == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Model for field " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return evaluator;
};


Teuchos::RCP<FieldEvaluator> State::GetFieldEvaluator_(Key key) {
  FieldEvaluatorMap::iterator lb = field_evaluators_.lower_bound(key);
  if (lb != field_evaluators_.end() && !(field_evaluators_.key_comp()(key, lb->first))) {
    return lb->second;
  } else {
    return Teuchos::null;
  }
};


void State::SetFieldEvaluator(Key key, const Teuchos::RCP<FieldEvaluator>& evaluator) {
  ASSERT(field_evaluators_[key] == Teuchos::null);
  field_evaluators_[key] = evaluator;
};


// -----------------------------------------------------------------------------
// State handles data management.
// -----------------------------------------------------------------------------
Teuchos::RCP<Field> State::GetField_(Key fieldname) {
  FieldMap::iterator lb = fields_.lower_bound(fieldname);
  if (lb != fields_.end() && !(fields_.key_comp()(fieldname, lb->first))) {
    return lb->second;
  } else {
    return Teuchos::null;
  }
};

Teuchos::RCP<const Field> State::GetField_(Key fieldname) const {
  FieldMap::const_iterator lb = fields_.lower_bound(fieldname);
  if (lb != fields_.end() && !(fields_.key_comp()(fieldname, lb->first))) {
    return lb->second;
  } else {
    return Teuchos::null;
  }
};


// Scalar requires
void State::RequireScalar(Key fieldname, Key owner) {
  Teuchos::RCP<Field> field = CheckConsistent_or_die_(fieldname, CONSTANT_SCALAR, owner);

  if (field == Teuchos::null) {
    Teuchos::RCP<FieldScalar> newfield = Teuchos::rcp(new FieldScalar(fieldname, owner));
    fields_[fieldname] = newfield;
  } else if (owner != Key("state")) {
    field->set_owner(owner);
  }
};


// Require a constant vector of a given size.
void State::RequireConstantVector(Key fieldname, Key owner,
              int dimension) {
  Teuchos::RCP<Field> field = CheckConsistent_or_die_(fieldname, CONSTANT_VECTOR, owner);

  if (field == Teuchos::null) {
    Teuchos::RCP<FieldConstantVector> field =
        Teuchos::rcp(new FieldConstantVector(fieldname, Key("state"), dimension));
    fields_[fieldname] = field;
  } else {
    Teuchos::RCP<FieldConstantVector> cv =
        Teuchos::rcp_static_cast<FieldConstantVector>(field);
    cv->set_dimension(dimension);

    if (owner != Key("state")) {
      cv->set_owner(owner);
    }
  }
};

void State::RequireConstantVector(Key fieldname, int dimension) {
  RequireConstantVector(fieldname, Key("state"), dimension);
};

// Vector Field requires
Teuchos::RCP<CompositeVectorFactory>
State::RequireField(Key fieldname, Key owner) {
  Teuchos::RCP<Field> field = CheckConsistent_or_die_(fieldname,
          COMPOSITE_VECTOR_FIELD, owner);

  if (field == Teuchos::null) {
    // Create the field and CV factory.
    Teuchos::RCP<FieldCompositeVector> field =
        Teuchos::rcp(new FieldCompositeVector(fieldname, owner));
    fields_[fieldname] = field;
    field_factories_[fieldname] = Teuchos::rcp(new CompositeVectorFactory());
  } else if (owner != Key("state")) {
    field->set_owner(owner);
  }

  return field_factories_[fieldname];
};


void State::RequireGravity() {
  int dim = 3;
  Key fieldname("gravity");
  RequireConstantVector(fieldname, dim);
  std::vector<Key> subfield_names(dim);
  subfield_names[0] = "x";
  if (dim > 1) subfield_names[1] = "y";
  if (dim > 2) subfield_names[2] = "z";
  Teuchos::RCP<Field> field = GetField_(fieldname);
  Teuchos::RCP<FieldConstantVector> cvfield =
    Teuchos::rcp_dynamic_cast<FieldConstantVector>(field, true);
  cvfield->set_subfield_names(subfield_names);
};


// -- access methods -- Const methods should be used by PKs who don't own
// the field, i.e.  flow accessing a temperature field if an energy PK owns
// the temperature field.  This ensures a PK cannot mistakenly alter data it
// doesn't own.  Non-const methods get used by the owning PK.
Teuchos::RCP<const double>
State::GetScalarData(Key fieldname) const {
  return GetField(fieldname)->GetScalarData();
};

Teuchos::RCP<double>
State::GetScalarData(Key fieldname, Key pk_name) {
  return GetField(fieldname, pk_name)->GetScalarData();
};

Teuchos::RCP<const Epetra_Vector>
State::GetConstantVectorData(Key fieldname) const {
  return GetField(fieldname)->GetConstantVectorData();
};

Teuchos::RCP<Epetra_Vector>
State::GetConstantVectorData(Key fieldname, Key pk_name) {
  return GetField(fieldname, pk_name)->GetConstantVectorData();
};

Teuchos::RCP<const CompositeVector>
State::GetFieldData(Key fieldname) const {
  return GetField(fieldname)->GetFieldData();
};

Teuchos::RCP<CompositeVector>
State::GetFieldData(Key fieldname, Key pk_name) {
  return GetField(fieldname, pk_name)->GetFieldData();
};


// Access to the full field record, not just the data.
Teuchos::RCP<Field> State::GetField(Key fieldname, Key pk_name) {
  Teuchos::RCP<Field> record = GetField_(fieldname);

  if (record == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Field " << fieldname << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  } else if (record->owner() != pk_name) {
    std::stringstream messagestream;
    messagestream << "PK " << pk_name
                  << " is attempting write access to field " << fieldname
                  << " which is owned by " << GetField_(fieldname)->owner();
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return record;
};

Teuchos::RCP<const Field> State::GetField(Key fieldname) const {
  Teuchos::RCP<const Field> record = GetField_(fieldname);

  if (record == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Field " << fieldname << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return record;
};

// modify methods
// -- modify by pointer, no copy
void State::SetData(Key fieldname, Key pk_name,
                    const Teuchos::RCP<double>& data) {
  GetField(fieldname, pk_name)->SetData(data);
};

void State::SetData(Key fieldname, Key pk_name,
                    const Teuchos::RCP<Epetra_Vector>& data) {
  GetField(fieldname, pk_name)->SetData(data);
};

void State::SetData(Key fieldname, Key pk_name,
                    const Teuchos::RCP<CompositeVector>& data){
  GetField(fieldname, pk_name)->SetData(data);
};


// -----------------------------------------------------------------------------
// State handles model parameters.
// -----------------------------------------------------------------------------
Teuchos::ParameterList
State::GetModelParameters(std::string modelname) {
  ASSERT(state_plist_.isSublist("model parameters"));
  Teuchos::ParameterList model_plist = state_plist_.sublist("model parameters");
  ASSERT(model_plist.isSublist(modelname));
  return model_plist.sublist(modelname);
}


// -----------------------------------------------------------------------------
// Initialization, etc.
// -----------------------------------------------------------------------------
// Initialize data, allowing values to be specified here or in the owning PK.
// All independent variables must be initialized here.
void State::Setup() {
  // State-required data
  if (state_plist_.isParameter("visualize mesh ranks")) {
    for (mesh_iterator mesh_it=meshes_.begin();
         mesh_it!=meshes_.end(); ++mesh_it) {
      Key rank_key;
      if (mesh_it->first == "domain") {
        rank_key = "mpi_comm_rank";
      } else {
        rank_key = mesh_it->first + "_mpi_comm_rank";
      }
      Teuchos::ParameterList plist;
      plist.set("evaluator name", rank_key);
      plist.set("mesh name", mesh_it->first);
      SetFieldEvaluator(rank_key, Teuchos::rcp(new RankModel(plist)));
    }
  }

  // Ensure compatibility of all the evaluators -- each evaluator's dependencies must
  // provide what is required of that evaluator.
  for (FieldEvaluatorMap::iterator evaluator=field_evaluators_.begin();
       evaluator!=field_evaluators_.end(); ++evaluator) {
    if (!evaluator->second->ProvidesKey(evaluator->first)) {
      std::stringstream messagestream;
      messagestream << "Field Evaluator " << evaluator->first
                    << " does not provide it's own key.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
    evaluator->second->EnsureCompatibility(Teuchos::ptr(this));
  }

  // Create all data for vector fields.
  // -- First use factories to instantiate composite vectors.
  for (FieldFactoryMap::iterator fac_it=field_factories_.begin();
       fac_it!=field_factories_.end(); ++fac_it) {
    GetField_(fac_it->first)->SetData(fac_it->second->CreateVector(true));
  }

  // -- Now create the data for all fields.
  for (field_iterator f_it = field_begin();
       f_it != field_end(); ++f_it) {
    f_it->second->CreateData();
  }
};


void State::Initialize() {
  // Initialize any other fields from state plist.
  InitializeFields_();

  // Ensure that non-evaluator-based fields are initialized.
  CheckNotEvaluatedFieldsInitialized_();

  // Initialize other field evaluators.
  InitializeEvaluators_();

  // Ensure everything is owned and initialized.
  CheckAllFieldsInitialized_();
};


void State::InitializeEvaluators_() {
  for (evaluator_iterator f_it = field_evaluator_begin();
       f_it != field_evaluator_end(); ++f_it) {
    f_it->second->HasFieldChanged(Teuchos::Ptr<State>(this), "state");
    fields_[f_it->first]->set_initialized();
  }
};


void State::InitializeFields_() {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    if (!f_it->second->initialized()) {
      if (state_plist_.isSublist("initial conditions")) {
        if (state_plist_.sublist("initial conditions").isSublist(f_it->first)) {
          Teuchos::ParameterList sublist = state_plist_.sublist("initial conditions").sublist(f_it->first);
          f_it->second->Initialize(sublist);
        }
      }
    }
  }
};


// Make sure all fields that are not evaluated by a FieldEvaluator are
// initialized.  Such fields may be used by an evaluator field but are not in
// the dependency tree due to poor design.
bool State::CheckNotEvaluatedFieldsInitialized_() {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    if (!HasFieldEvaluator(f_it->first) && !field->initialized()) {
      std::stringstream messagestream;
      messagestream << "Field " << field->fieldname() << " was not initialized.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
      return false;
    }
  }
  return true;
};


// Make sure all fields have gotten their IC, either from State or the owning PK.
bool State::CheckAllFieldsInitialized_() {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    if (!field->initialized()) {
      // field was not initialized
      std::stringstream messagestream;
      messagestream << "Field " << field->fieldname() << " was not initialized.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
      return false;
    }
  }
  return true;
};


// Check the consistency of a field's meta-data.
Teuchos::RCP<Field>
State::CheckConsistent_or_die_(Key fieldname, FieldType type, Key owner) {
  ASSERT(fieldname != Key(""));

  Teuchos::RCP<Field> record = GetField_(fieldname);
  if (record == Teuchos::null) {
    // No existing field of that fieldname, so anything is ok.
    return record;
  }

  if (owner == Key("state") ||
      record->owner() == Key("state") ||
      owner == record->owner()) {

    if (record->type() != type) {
      std::stringstream messagestream;
      messagestream << "Requested field " << fieldname << " of type "
                    << type << " already exists and is of type "
                    << record->type();
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
  } else {
    // Field exists and is owned
    std::stringstream messagestream;
    messagestream << "Requested field " << fieldname << " already exists and is owned by "
                  << record->owner();
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  return record;
};


// Set the State's time and evaluate any state-owned scalar/vector fields at
// the new time if they have an associated function.
void State::set_time( double new_time ) {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    f_it->second->Compute(new_time);
  }
  time_ = new_time;
}


// Non-member function for vis.
void WriteVis(const Teuchos::Ptr<Visualization>& vis,
              const Teuchos::Ptr<State>& S) {
  if (!vis->is_disabled()) {
    // Create the new time step (internally we use seconds, but write the time in years).
    vis->WriteMesh(S->time()/(365.25*24*60*60),S->cycle());
    vis->CreateTimestep(S->time()/(365.25*24*60*60),S->cycle());

    // Write all fields to the visualization file, the fields know if they
    // need to be written.
    for (State::field_iterator field=S->field_begin(); field!=S->field_end(); ++field) {
      field->second->WriteVis(vis);
    }

    // Finalize i/o.
    vis->FinalizeTimestep();
  }
};


// Non-member function for checkpointing.
void WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& chk,
                     const Teuchos::Ptr<State>& S,
                     double dt) {
  if ( !chk->is_disabled() ) {
    chk->CreateFile(S->cycle());

    for (State::field_iterator field=S->field_begin(); field!=S->field_end(); ++field) {
      field->second->WriteCheckpoint(chk);
    }

    chk->WriteAttributes(S->time(), dt, S->cycle());
  }
};

// Non-member function for checkpointing.
double ReadCheckpoint(Epetra_MpiComm* comm,
                      const Teuchos::Ptr<State>& S,
                      std::string filename) {
  Teuchos::Ptr<HDF5_MPI> checkpoint = Teuchos::ptr(new HDF5_MPI(*comm, filename));

  // load the attributes
  double time(0.);
  checkpoint->readAttrReal(time, "time");
  S->set_time(time);

  double dt(0.);
  checkpoint->readAttrReal(dt, "dt");

  int cycle(0);
  checkpoint->readAttrInt(cycle, "cycle");
  S->set_cycle(cycle);

  // load the number of processes and ensure they are the same -- otherwise
  // the below just gives crap.
  int rank(-1);
  checkpoint->readAttrInt(rank, "mpi_comm_world_rank");
  if (comm->NumProc() != rank) {
    std::stringstream messagestream;
    messagestream << "Requested checkpoint file " << filename << " was created on "
                  << rank << " processes, making it incompatible with this run on "
                  << comm->NumProc() << " process.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // load the data
  for (State::field_iterator field=S->field_begin(); field!=S->field_end(); ++field) {
    if (field->second->io_checkpoint()) {
      field->second->ReadCheckpoint(checkpoint);
    }
  }

  return dt;
};

// Non-member function for checkpointing.
double ReadCheckpointInitialTime(Epetra_MpiComm* comm,
        std::string filename) {
  Teuchos::Ptr<HDF5_MPI> checkpoint = Teuchos::ptr(new HDF5_MPI(*comm, filename));

  // load the attributes
  double time(0.);
  checkpoint->readAttrReal(time, "time");
  return time;
};

// Non-member function for deforming the mesh after reading a checkpoint file
// that contains the vertex coordinate field (this is written by deformation pks)
void DeformCheckpointMesh(const Teuchos::Ptr<State>& S) {
  if (S->HasField("vertex coordinate")) { // only deform mesh if vertex coordinate field exists
    AmanziMesh::Mesh * write_access_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S->GetMesh());
    // get vertex coordinates state
    const Epetra_MultiVector& vc = *S->GetFieldData("vertex coordinate")
        ->ViewComponent("node",false);
    int dim = write_access_mesh_->space_dimension();
    int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
                                              Amanzi::AmanziMesh::OWNED);  
    Amanzi::AmanziMesh::Entity_ID_List nodeids;  
    Amanzi::AmanziGeometry::Point new_coords(3);
    AmanziGeometry::Point_List new_pos, final_pos;  
    // loop over vertices and update vc
    for (int iV=0; iV<nV; iV++) {
      // set the coords of the node
      for ( int s=0; s<dim; ++s ) {
        new_coords[s] = vc[s][iV];
      }
      // push back for deform method
      nodeids.push_back(iV);
      new_pos.push_back(new_coords);    
    } 
    write_access_mesh_->deform( nodeids, new_pos, true, &final_pos); // deforms the mesh
  }
  
}

} // namespace amanzi
