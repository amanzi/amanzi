/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! State is the primary data manager for Amanzi-ATS

/*!
  Interface for the State.  State is a simple data-manager, allowing PKs to
  require, read, and write various fields, including:

  -- Acts as a factory for fields through the various require methods.
  -- Provides some data protection by providing both const and non-const
       data pointers to PKs.

*/

#include <iostream>
#include <map>
#include <ostream>

#include "Epetra_Vector.h"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "CompositeVector.hh"
#include "errors.hh"

#include "State.hh"
#include "evaluator/Evaluator_Factory.hh"

namespace Amanzi {

State::State() : state_plist_(Teuchos::rcp(new Teuchos::ParameterList())){};

State::State(const Teuchos::RCP<Teuchos::ParameterList>& state_plist)
  : state_plist_(state_plist){};

// -----------------------------------------------------------------------------
// State handles mesh management.
// -----------------------------------------------------------------------------
void
State::RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                          bool deformable)
{
  RegisterMesh("domain", mesh, deformable);
}

void
State::RegisterMesh(Key key, Teuchos::RCP<AmanziMesh::Mesh> mesh,
                    bool deformable)
{
  meshes_.emplace(std::piecewise_construct,
                  std::forward_as_tuple<Key>(std::move(key)),
                  std::forward_as_tuple<Teuchos::RCP<AmanziMesh::Mesh>, bool>(
                    std::move(mesh), std::move(deformable)));
};

void
State::AliasMesh(const Key& target, Key alias)
{
  bool deformable = IsDeformableMesh(target);

  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  try {
    mesh = meshes_.at(target).first;
  } catch (const std::out_of_range& e) {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << target << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }

  RegisterMesh(alias, mesh, deformable);
};

void
State::RemoveMesh(const Key& key)
{
  meshes_.erase(key);
};

bool
State::HasMesh(const Key& key) const
{
  return Keys::hasKey(meshes_, key);
}

Teuchos::RCP<const AmanziMesh::Mesh>
State::GetMesh(const Key& key) const
{
  try {
    return meshes_.at(key).first;
  } catch (const std::out_of_range& e) {
    std::stringstream messagestream;
    messagestream << "Mesh \"" << key << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
};

Teuchos::RCP<AmanziMesh::Mesh>
State::GetDeformableMesh(const Key& key)
{
  std::pair<Teuchos::RCP<AmanziMesh::Mesh>, bool> m;
  try {
    m = meshes_.at(key);
  } catch (const std::out_of_range& e) {
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

bool
State::IsDeformableMesh(const Key& key) const
{
  try {
    return meshes_.at(key).second;
  } catch (const std::out_of_range& e) {
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
Evaluator&
State::RequireEvaluator(const Key& key, const Key& tag)
{
  // does it already exist?
  if (HasEvaluator(key, tag)) return GetEvaluator(key, tag);

  // See if the key is provided by another existing evaluator.
  for (auto& f_it : evaluators_) {
    if (Keys::hasKey(f_it.second, tag) &&
        f_it.second.at(tag)->ProvidesKey(key, tag)) {
      auto& evaluator = f_it.second.at(tag);
      SetEvaluator(key, tag, evaluator);
      return *evaluator;
    }
  }

  // Create the evaluator from State's plist
  // -- Get the Field Evaluator plist
  Teuchos::ParameterList& fm_plist = state_plist_->sublist("evaluators");

  if (fm_plist.isSublist(key)) {
    // -- Get _a_copy_ of this evaluator's plist
    Teuchos::ParameterList sublist(fm_plist.sublist(key));
    sublist.set("tag", tag);

    // -- Insert any model parameters.
    if (sublist.isParameter("common model parameters")) {
      std::string modelname = sublist.get<std::string>("common model parameters");
      Teuchos::ParameterList modellist = GetModelParameters(modelname);
      sublist.set("model parameters", modellist);
    // } else if (sublist.isParameter("models parameters")) {
    //   Teuchos::Array<std::string> modelnames =
    //     sublist.get<Teuchos::Array<std::string>>("models parameters");
    //   for (Teuchos::Array<std::string>::const_iterator modelname =
    //          modelnames.begin();
    //        modelname != modelnames.end();
    //        ++modelname) {
    //     Teuchos::ParameterList modellist = GetModelParameters(*modelname);
    //     std::string modeltype = modellist.get<std::string>("model type");
    //     sublist.set(modeltype, modellist);
    //   }
    }

    // -- Create and set the evaluator.
    Evaluator_Factory evaluator_factory;
    auto evaluator = evaluator_factory.createEvaluator(sublist);
    SetEvaluator(key, tag, evaluator);
    return *evaluator;
  }

  // cannot find the evaluator, error
  Errors::Message message;
  message << "Model for field " << key << " cannot be created in State.";
  throw(message);
}

bool
State::HasEvaluator(const Key& key, const Key& tag)
{
  if (Keys::hasKey(evaluators_, key)) {
    return Keys::hasKey(evaluators_.at(key), tag);
  } else {
    return false;
  }
}

Evaluator&
State::GetEvaluator(const Key& key, const Key& tag)
{
  return *GetEvaluatorPtr(key, tag);
}

const Evaluator&
State::GetEvaluator(const Key& key, const Key& tag) const
{
  try {
    return *evaluators_.at(key).at(tag);
  } catch (std::out_of_range) {
    std::stringstream messagestream;
    messagestream << "Evaluator for field \"" << key << "\" at tag \"" << tag
                  << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
};

Teuchos::RCP<Evaluator>
State::GetEvaluatorPtr(const Key& key, const Key& tag)
{
  try {
    return evaluators_.at(key).at(tag);
  } catch (std::out_of_range) {
    std::stringstream messagestream;
    messagestream << "Evaluator for field \"" << key << "\" at tag \"" << tag
                  << "\" does not exist in the state.";
    Errors::Message message(messagestream.str());
    throw(message);
  }
};

void
State::SetEvaluator(const Key& key, const Key& tag,
                    const Teuchos::RCP<Evaluator>& evaluator)
{
  evaluators_[key][tag] = evaluator;
  for (const auto& dep_tag : evaluator->dependencies()) {
    RequireEvaluator(dep_tag.first, dep_tag.second);
  }  
};

void
State::WriteDependencyGraph() const
{
  std::ofstream os("dependency_graph.txt", std::ios::out);
  for (auto& fe : evaluators_) {
    for (auto& e : fe.second) { os << *e.second; }
  }
  os.close();
}

// -----------------------------------------------------------------------------
// State handles model parameters.
// -----------------------------------------------------------------------------
Teuchos::ParameterList
State::GetModelParameters(std::string modelname)
{
  if (!state_plist_->isSublist("common model parameters")) {
    Errors::Message message(
      "State parameter list is missing \"common model parameters\" sublist.");
    throw(message);
  }
  Teuchos::ParameterList& model_plist =
    state_plist_->sublist("common model parameters");
  if (!model_plist.isSublist(modelname)) {
    Errors::Message message;
    message << "State \"model parameters\" list does not have a model named \""
            << modelname << "\".";
    throw(message);
  }
  return model_plist.sublist(modelname);
}

// -----------------------------------------------------------------------------
// Initialization, etc.
// -----------------------------------------------------------------------------
// Initialize data, allowing values to be specified here or in the owning PK.
// All independent variables must be initialized here.
void
State::Setup()
{
  Require<double>("time", "", "time");
  Require<double>("time", "next", "time");
  Require<int>("cycle", "", "cycle");
  Require<int>("cycle", "next", "cycle");

  // Ensure compatibility of all the evaluators -- each evaluator's dependencies
  // must provide what is required of that evaluator.
  for (auto& e : evaluators_) {
    for (auto& r : e.second) {
      if (!r.second->ProvidesKey(e.first, r.first)) {
        AMANZI_ASSERT(r.second->ProvidesKey(e.first, r.first));
        Errors::Message message;
        message << "Evaluator \"" << e.first << "\" with tag \"" << r.first
                << "\" does not provide its own key.";
        throw(message);
      }
      r.second->EnsureCompatibility(*this);
    }
  }

  // // Now make sure derivatives are properly set up as well
  // for (const auto& deriv : derivs_) {
  //   auto keytag = Keys::splitKeyTag(deriv.first);
  //   for (const auto& wrt : *deriv.second) {
  //     auto wrtkeytag = Keys::splitKeyTag(wrt.first);
  //     GetEvaluator(keytag.first,keytag.second)
  //         .EnsureCompatibleDerivative(*this,wrtkeytag.first,wrtkeytag.second);
  //   }
  // }

  // -- Create the data for all fields.
  for (auto& f : data_) { f.second->CreateData(); }
  for (auto& deriv : derivs_) { deriv.second->CreateData(); }

  Set("time", "", "time", 0.0);
  GetRecordW("time", "time").set_initialized();
  Set("time", "next", "time", 0.0);
  GetRecordW("time", "next", "time").set_initialized();
  Set("cycle", "", "cycle", 0);
  GetRecordW("cycle", "cycle").set_initialized();
  Set("cycle", "next", "cycle", 0);
  GetRecordW("cycle", "next", "cycle").set_initialized();
};

void
State::Initialize()
{
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
void
State::Initialize(const State& other)
{
  // initialize data in this from data in other
  for (auto& e : data_) {
    if (!e.second->initialized()) {
      for (auto& et : *e.second) {
        auto& fieldname = e.first;
        auto& tag = et.first;
        auto& record = *et.second;
        if (other.HasData(fieldname, tag)) {
          record = other.GetRecord(fieldname, tag);
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

void
State::InitializeEvaluators(){
  // for (auto& f_it : evaluators_) {
  //   for (auto& e : f_it.second) {
  //     std::cout << "Initing eval: \"" << f_it.first << "\" tag \"" <<
  //     e.first << "\"" << std::endl;
  //     e.second->Update(*this, "state");
  //     GetRecordW(f_it.first, e.first, f_it.first).set_initialized();
  //   }
  // }
};

void
State::InitializeFields()
{
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
void
State::AliasEvaluators()
{
  for (auto& f_it : data_) {
    auto fname = f_it.first;
    for (auto& e : *f_it.second) {
      auto tagname = e.first;
      if (!HasEvaluator(fname, tagname)) {
        // first check and see if there is a Evaluator, but we haven't yet used
        // it.
        Teuchos::RCP<Evaluator> found_eval;
        for (auto& f_eval_it : evaluators_) {
          if (Keys::hasKey(f_eval_it.second, tagname) &&
              f_eval_it.second.at(tagname)->ProvidesKey(fname, tagname)) {
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
bool
State::CheckNotEvaluatedFieldsInitialized()
{
  // for (auto& f_it : data_) {
  //   for (auto& field : *f_it.second) {
  //     if (!HasEvaluator(f_it.first, field.first) &&
  //     !field.second->initialized()) {
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

// Make sure all fields have gotten their IC, either from State or the owning
// PK.
bool
State::CheckAllFieldsInitialized()
{
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
void
WriteVis(Visualization& vis, const State& S)
{
  if (!vis.is_disabled()) {
    // Create the new time step
    // NOTE: internally we use seconds, but for vis we write the time in years
    vis.CreateFile(S.time(), S.cycle());

    // Write all fields to the visualization file, the fields know if they
    // need to be written.
    for (auto r = S.data_begin(); r != S.data_end(); ++r) {
      r->second->WriteVis(vis);
    }

    // vis.WriteRegions();
    // vis.WritePartition();

    // Finalize i/o.
    vis.FinalizeFile();
  }
};

// Non-member function for checkpointing.
void
WriteCheckpoint(Checkpoint& chkp, const State& S, bool final)
{
  chkp.CreateFile(S.time(), S.cycle());
  for (auto r = S.data_begin(); r != S.data_end(); ++r) {
    r->second->WriteCheckpoint(chkp);
  }
  chkp.FinalizeFile(final);
};

// Non-member function for checkpointing.
void
ReadCheckpoint(const Comm_ptr_type& comm, State& S, const std::string& filename)
{
  Teuchos::ParameterList plist;
  plist.set("file name", filename);
  plist.set("file type", "HDF5");
  Checkpoint chkp(plist, comm, true);

  // Load the number of processes and ensure they are the same.
  int num_procs(-1);

  Teuchos::ParameterList mpi_num_procs("mpi_num_procs");
  chkp.Read(mpi_num_procs, num_procs);
  if (comm->getSize() != num_procs) {
    std::stringstream messagestream;
    messagestream << "Requested checkpoint file " << filename
                  << " was created on " << num_procs
                  << " processes, making it incompatible with this run on "
                  << comm->getSize() << " processes.";
    Errors::Message message(messagestream.str());
    throw(message);
  }

  // load the data
  for (auto data = S.data_begin_(); data != S.data_end_(); ++data) {
    data->second->ReadCheckpoint(chkp);
  }
};

// Non-member function for deforming the mesh after reading a checkpoint file
// that contains the vertex coordinate field (this is written by deformation
// pks)
//
// FIX ME: Refactor this to make the name more general.  Should align with a
// mesh name prefix or something, and the coordinates should be written by
// state in WriteCheckpoint if mesh IsDeformableMesh() --ETC
void
DeformCheckpointMesh(State& S)
{
  // FIXME EPETRA TO TPETRA: deforming
  // if (S.HasData("vertex coordinate")) { // only deform mesh if vertex
  // coordinate
  //                                       // field exists
  //   AmanziMesh::Mesh *write_access_mesh =
  //       const_cast<AmanziMesh::Mesh *>(&*S.GetMesh());

  //   // get vertex coordinates state
  //   const CompositeVector &vc = S.Get<CompositeVector>("vertex coordinate");
  //   vc.ScatterMasterToGhosted("node");
  //   auto vc_n = vc.ViewComponent("node", true);

  //   int dim = write_access_mesh->space_dimension();
  //   Amanzi::AmanziMesh::Entity_ID_List nodeids;
  //   Amanzi::AmanziGeometry::Point new_coords(dim);
  //   AmanziGeometry::Point_List new_pos, final_pos;
  //   // loop over vertices and update vc
  //   unsigned int nV = vc_n.getLocalLength();
  //   for (unsigned int iV = 0; iV != nV; ++iV) {
  //     // set the coords of the node
  //     for (unsigned int s = 0; s != dim; ++s)
  //       new_coords[s] = vc_n[s][iV];

  //     // push back for deform method
  //     nodeids.push_back(iV);
  //     new_pos.push_back(new_coords);
  //   }

  //   // deform
  //   write_access_mesh->deform(nodeids, new_pos, true,
  //                             &final_pos); // deforms the mesh
  //  }
}

} // namespace Amanzi
