/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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
#include <map>
#include <ostream>

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "CompositeVector.hh"

#include "State.hh"
#include "evaluator/Evaluator_Factory.hh"

namespace Amanzi {

State::State()
    : state_plist_(Teuchos::rcp(new Teuchos::ParameterList()))
{};

State::State(const Teuchos::RCP<Teuchos::ParameterList>& state_plist)
    : state_plist_(state_plist)
{};


// -----------------------------------------------------------------------------
// State handles mesh management.
// -----------------------------------------------------------------------------
void State::RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                               bool deformable) {
  RegisterMesh("domain", mesh, deformable);
}


void State::RegisterMesh(Key key,
                         Teuchos::RCP<AmanziMesh::Mesh> mesh,
                         bool deformable) {
    meshes_.emplace(std::piecewise_construct,
                    std::forward_as_tuple<Key>(std::move(key)),
                    std::forward_as_tuple<Teuchos::RCP<AmanziMesh::Mesh>, bool>(std::move(mesh),std::move(deformable)));
};


void State::AliasMesh(const Key& target, Key alias) {
  bool deformable = IsDeformableMesh(target);

  Teuchos::RCP<AmanziMesh::Mesh> mesh;  
  try {
    mesh = meshes_.at(target).first;
  } catch(const std::out_of_range& e) {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << target << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }

  RegisterMesh(alias, mesh, deformable);
};


void State::RemoveMesh(const Key& key) {
  meshes_.erase(key);
};


bool State::HasMesh(const Key& key) const {
  return Keys::hasKey(meshes_, key);
}

Teuchos::RCP<const AmanziMesh::Mesh> State::GetMesh(const Key& key) const {
  try {
    return meshes_.at(key).first;
  } catch(const std::out_of_range& e) {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << key << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
};


Teuchos::RCP<AmanziMesh::Mesh> State::GetDeformableMesh(const Key& key) {
  std::pair<Teuchos::RCP<AmanziMesh::Mesh>,bool> m;
  try {
    m = meshes_.at(key);
  } catch(const std::out_of_range& e) {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << key << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }

  if (m.second) {
    return m.first;
  } else {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << key << "\" is not deformable.";
    Errors::Message message(messagestream.str());
    throw(message);
  }

};


bool State::IsDeformableMesh(const Key& key) const {
  try {
    return meshes_.at(key).second;
  } catch(const std::out_of_range& e) {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << key << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
  return false;
};


// -----------------------------------------------------------------------------
// State handles data evaluation.
// -----------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
State::RequireEvaluator(const Key& key, const Key& tag) {

  // does it already exist?
  if (HasEvaluator(key, tag)) {
    return GetEvaluator(key, tag);
  }

  // See if the key is provided by another existing evaluator.
  for (auto& f_it : evaluators_) {
    if (Keys::hasKey(f_it.second, tag) &&
        f_it.second.at(tag)->ProvidesKey(key)) {
      auto& evaluator = f_it.second.at(tag);
      SetEvaluator(key, evaluator);
      return evaluator;
    }
  }

  // Create the evaluator from State's plist
  // -- Get the Field Evaluator plist
  Teuchos::ParameterList& fm_plist =
      state_plist_->sublist("evaluators");

  if (fm_plist.isSublist(key)) {
    // -- Get this evaluator's plist.
    Teuchos::ParameterList sublist = fm_plist.sublist(key);

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
    Evaluator_Factory evaluator_factory;
    sublist.set("tag", tag);
    auto evaluator = evaluator_factory.createEvaluator(sublist);
    SetEvaluator(key, tag, evaluator);
    return evaluator;
  }

  // cannot find the evaluator, error
  std::stringstream messagestream;
  messagestream << "Model for field " << key << " cannot be created in State.";
  Errors::Message message(messagestream.str());
  throw(message);
  return Teuchos::null;
}

bool State::HasEvaluator(const Key& key, const Key& tag) {
  if (Keys::hasKey(evaluators_,key) ) {
    return Keys::hasKey(evaluators_.at(key), tag);
  } else {
    return false;
  }  
}

Teuchos::RCP<Evaluator> State::GetEvaluator(const Key& key, const Key& tag) {
  try {
    return evaluators_.at(key).at(tag);
  } catch(std::out_of_range) {
    std::stringstream messagestream;
    messagestream << "Evaluator for field \"" << key << "\" at tag \"" << tag
                  << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
};


Teuchos::RCP<const Evaluator> State::GetEvaluator(const Key& key,
        const Key& tag) const {
  try {
    return evaluators_.at(key).at(tag);
  } catch(std::out_of_range) {
    std::stringstream messagestream;
    messagestream << "Evaluator for field \"" << key << "\" at tag \"" << tag
                  << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
};

void State::SetEvaluator(const Key& key, const Key& tag, const Teuchos::RCP<Evaluator>& evaluator) {
  evaluators_[key][tag] = evaluator;
};

void State::WriteDependencyGraph() const {
  std::ofstream os("dependency_graph.txt", std::ios::out);
  for (auto& fe : evaluators_) {
    for (auto& e : fe.second) {
       os << *e.second;
    }
  }
  os.close();
}


// -----------------------------------------------------------------------------
// State handles model parameters.
// -----------------------------------------------------------------------------
Teuchos::ParameterList
State::GetModelParameters(std::string modelname) {
  if (!state_plist_->isSublist("model parameters")) {
    Errors::Message message("State parameter list is missing \"model parameters\" sublist.");
    throw(message);
  }
  Teuchos::ParameterList& model_plist = state_plist_->sublist("model parameters");
  if (!model_plist.isSublist(modelname)) {
    Errors::Message message;
    message << "State \"model parameters\" list does not have a model named \"" << modelname << "\".";
    throw(message);
  }
  return model_plist.sublist(modelname);
}


// -----------------------------------------------------------------------------
// Initialization, etc.
// -----------------------------------------------------------------------------
// Initialize data, allowing values to be specified here or in the owning PK.
// All independent variables must be initialized here.
void State::Setup() {
  Require<double>("time", "", "time");
  Require<double>("time", "next", "time");
  Require<int>("cycle", "", "cycle");
  Require<int>("cycle", "next", "cycle");

  // Ensure compatibility of all the evaluators -- each evaluator's dependencies must
  // provide what is required of that evaluator.
  for (auto& e : evaluators_) {
    for (auto& r : e.second) {
      if (!r.second->ProvidesKey(e.first)) {
        std::stringstream messagestream;
        messagestream << "Evaluator \"" << e.first << "\" with tag \""
                      << r.first << "\" does not provide its own key.";
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
      }
      r.second->EnsureCompatibility(*this);
    }
  }

  // -- Create the data for all fields.
  for (auto& f : data_) {
    f.second->CreateData();
  }

  Set("time", "", "time", 0.0);
  GetRecordW("time", "time").set_initialized();
  Set("time", "next", "time", 0.0);
  GetRecordW("time", "next", "time").set_initialized();
  Set("cycle", "", "cycle", 0);
  GetRecordW("cycle", "cycle").set_initialized();
  Set("cycle", "next", "cycle", 0);
  GetRecordW("cycle", "next", "cycle").set_initialized();
};


void State::Initialize() {
  // Initialize any other fields from state plist.
  InitializeFields();

  // Alias evaluators: multiple fields may share an evaluator
  AliasEvaluators();
  
  // Ensure that non-evaluator-based fields are initialized.
  CheckNotEvaluatedFieldsInitialized();

  // Initialize other field evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  CheckAllFieldsInitialized();

  // Write dependency graph.
  WriteDependencyGraph();
};

// Initialized this from other state
void State::Initialize(const State& other) {

  // initialize data in this from data in other
  for (auto& e : data_) {
    if (!e.second->initialized()) {
      for (auto& et : *e.second) {
        auto& fieldname = e.first;
        auto& tag = et.first;
        auto& record = *et.second;
        if (other.HasData(fieldname, tag)) {
          record = other.GetRecord(fieldname,tag);
          record.set_initialized(true);
        }
      }
    }
  }

  // Initialize any other fields from state plist.
  InitializeFields();

  // Alias evaluators: multiple fields may share an evaluator
  AliasEvaluators();
  
  // Ensure that non-evaluator-based fields are initialized.
  CheckNotEvaluatedFieldsInitialized();

  // Initialize other field evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  CheckAllFieldsInitialized();

  // Write dependency graph.
  //
  // FIXME: This needs a stage identifyier or something?  Currently just
  // overwrites the last State's WriteDependencyGraph() output. --etc
  WriteDependencyGraph();
};

void State::InitializeEvaluators() {
  // for (auto& f_it : evaluators_) {
  //   for (auto& e : f_it.second) {
  //     std::cout << "Initing eval: \"" << f_it.first << "\" tag \"" << e.first << "\"" << std::endl;
  //     e.second->Update(*this, "state");
  //     GetRecordW(f_it.first, e.first, f_it.first).set_initialized();
  //   }
  // }
};


void State::InitializeFields() {
  if (state_plist_->isSublist("initial conditions")) {
    for (auto& e : data_) {
      if (!e.second->initialized()) {
        if (state_plist_->sublist("initial conditions").isSublist(e.first)) {
          Teuchos::ParameterList sublist =
              state_plist_->sublist("initial conditions").sublist(e.first);
          e.second->Initialize(sublist);
        }
      }
    }
  }
}


// Some fields share evaluators, these must be aliased by pointer.
void State::AliasEvaluators() {
  for (auto& f_it : data_) {
    auto fname = f_it.first;
    for (auto& e : *f_it.second) {
      auto tagname = e.first;
      if (!HasEvaluator(fname, tagname)) {
        // first check and see if there is a Evaluator, but we haven't yet used it.
        Teuchos::RCP<Evaluator> found_eval;
        for (auto& f_eval_it : evaluators_) {
          if (Keys::hasKey(f_eval_it.second, tagname) &&
              f_eval_it.second.at(tagname)->ProvidesKey(fname)) {
            found_eval = f_eval_it.second.at(tagname);
            break;
          }
        }

        if (found_eval.get()) {
          // found an evaluator that provides this key, set it.
          SetEvaluator(fname, tagname, found_eval);
        }
      }
    }
  }
}


// Make sure all fields that are not evaluated by a Evaluator are
// initialized.  Such fields may be used by an evaluator field but are not in
// the dependency tree due to poor design.
bool State::CheckNotEvaluatedFieldsInitialized() {
  // for (auto& f_it : data_) {
  //   for (auto& field : *f_it.second) {
  //     if (!HasEvaluator(f_it.first, field.first) && !field.second->initialized()) {
  //       // No evaluator, not intialized... FAIL.
  //       Errors::Message message;
  //       message << "Field \"" << f_it.first << "\" at tag \"" << field.first
  //               << "\" was not initialized.";
  //       throw(message);
  //       return false;
  //     }
  //   }
  // }
  return true;
};


// Make sure all fields have gotten their IC, either from State or the owning PK.
bool State::CheckAllFieldsInitialized() {
  // for (auto& f_it : data_) {
  //   Record& field = f_it.second->GetRecord("");
  //   if (!field.initialized()) {
  //     // field was not initialized
  //     std::stringstream messagestream;
  //     messagestream << "Field \"" << f_it.first << "\" was not initialized.";
  //     Errors::Message message(messagestream.str());
  //     throw(message);
  //     return false;
  //   }
  // }
  return true;
};


// Non-member function for vis.
void WriteVis(Visualization& vis,
              const State& S)
{
  if (!vis.is_disabled()) {
    // Create the new time step
    // NOTE: internally we use seconds, but for vis we write the time in years
    vis.CreateTimestep(S.time()/(365.25*24*60*60),S.cycle());

    // Write all fields to the visualization file, the fields know if they
    // need to be written.
    for (auto r=S.data_begin(); r!=S.data_end(); ++r) {
      r->second->WriteVis(vis);
    }
    
    vis.WriteRegions();
    vis.WritePartition();

    // Finalize i/o.
    vis.FinalizeTimestep();
  }
};


// Non-member function for checkpointing.
void WriteCheckpoint(Checkpoint& chkp,
                     const Epetra_MpiComm& comm,
                     const State& S,
                     bool final)
{
  if (!chkp.is_disabled()) {
    chkp.SetFinal(final);
    chkp.CreateFile(S.cycle());
    for (auto r=S.data_begin(); r!=S.data_end(); ++r) {
      r->second->WriteCheckpoint(chkp);
    }
    chkp.Write("mpi_num_procs", comm.NumProc());
    chkp.Finalize();
  }
};


// Non-member function for checkpointing.
void ReadCheckpoint(const Epetra_MpiComm& comm,
                    State& S,
                    const std::string& filename)
{
  Checkpoint chkp(filename, comm);

  // Load the number of processes and ensure they are the same.
  int num_procs(-1);
  chkp.Read("mpi_num_procs", num_procs);
  if (comm.NumProc() != num_procs) {
    std::stringstream messagestream;
    messagestream << "Requested checkpoint file " << filename << " was created on "
                  << num_procs << " processes, making it incompatible with this run on "
                  << comm.NumProc() << " processes.";
    Errors::Message message(messagestream.str());
    throw(message);
  }

  // load the data
  for (auto data=S.data_begin_(); data!=S.data_end_(); ++data) {
    data->second->ReadCheckpoint(chkp);
  }
  
  chkp.Finalize();
};


// Non-member function for deforming the mesh after reading a checkpoint file
// that contains the vertex coordinate field (this is written by deformation pks)
//
// FIX ME: Refactor this to make the name more general.  Should align with a
// mesh name prefix or something, and the coordinates should be written by
// state in WriteCheckpoint if mesh IsDeformableMesh() --ETC
void DeformCheckpointMesh(State& S) {
  if (S.HasData("vertex coordinate")) { // only deform mesh if vertex coordinate field exists
    AmanziMesh::Mesh* write_access_mesh =
        const_cast<AmanziMesh::Mesh*>(&*S.GetMesh());

    // get vertex coordinates state
    const CompositeVector& vc =
      S.Get<CompositeVector>("vertex coordinate");
    vc.ScatterMasterToGhosted("node");
    const Epetra_MultiVector& vc_n = *vc.ViewComponent("node",true);

    int dim = write_access_mesh->space_dimension();
    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    Amanzi::AmanziGeometry::Point new_coords(dim);
    AmanziGeometry::Point_List new_pos, final_pos;
    // loop over vertices and update vc
    unsigned int nV = vc_n.MyLength();
    for (unsigned int iV=0; iV!=nV; ++iV) {
      // set the coords of the node
      for (unsigned int s=0; s!=dim; ++s) new_coords[s] = vc_n[s][iV];

      // push back for deform method
      nodeids.push_back(iV);
      new_pos.push_back(new_coords);
    }

    // deform
    write_access_mesh->deform( nodeids, new_pos, true, &final_pos); // deforms the mesh
  }
}


} // namespace amanzi
