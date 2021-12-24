/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Implementation for the State.  State is a simple data-manager, allowing PKs to
  require, read, and write various fields.  Provides some data protection by
  providing both const and non-const fields to PKs.  Provides some
  initialization capability -- this is where all independent variables can be
  initialized (as independent variables are owned by state, not by any PK).
*/

#include <iostream>
#include <map>
#include <ostream>
#include <regex>

#include "boost/algorithm/string/predicate.hpp"
#include "boost/filesystem.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_Vector.h"

// Amanzi
#include "CompositeVector.hh"
#include "DomainSet.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "MeshPartition.hh"
#include "StringExt.hh"

// Amanzi::State
#include "Evaluator_Factory.hh"
#include "EvaluatorCellVolume.hh"
#include "EvaluatorPrimary.hh"
#include "State.hh"
#include "StateDefs.hh"

namespace Amanzi {

State::State() {
  vo_ = Teuchos::rcp(new VerboseObject("State", "none"));
};

State::State(Teuchos::ParameterList& state_plist) :
    state_plist_(state_plist),
    time_(0.0),
    cycle_(0),
    position_in_tp_(TIME_PERIOD_START) {
  vo_ = Teuchos::rcp(new VerboseObject("State", state_plist_));
};


// -----------------------------------------------------------------------------
// State handles mesh management.
// -----------------------------------------------------------------------------
void State::RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                               bool deformable) {
  RegisterMesh("domain", mesh, deformable);
}


void State::RegisterMesh(const Key& key,
                         const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                         bool deformable) {
  meshes_.insert(std::make_pair(key, std::make_pair(mesh,deformable)));
};


Teuchos::RCP<const AmanziMesh::Mesh> State::GetMesh(const Key& key) const
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  if (key.empty()) {
    mesh = GetMesh_("domain");
  } else {
    mesh = GetMesh_(key);
  }
  if (mesh == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Mesh " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return mesh;
};


Teuchos::RCP<AmanziMesh::Mesh> State::GetDeformableMesh(Key key)
{
  if (key.empty()) key = "domain";

  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    if (lb->second.second) {
      return lb->second.first;
    } else {
      std::stringstream messagestream;
      messagestream << "Mesh " << key << " is not deformable.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
  } else {
    std::stringstream messagestream;
    messagestream << "Mesh " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return Teuchos::null;
}


bool State::IsDeformableMesh(const Key& key) const
{
  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second.second;
  } else {
    Errors::Message msg;
    msg << "Mesh " << key << " does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  return false;
}


Teuchos::RCP<AmanziMesh::Mesh> State::GetMesh_(const Key& key) const
{
  if (key.empty()) return GetMesh_("domain");
  
  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second.first;
  } else {
    return Teuchos::null;
  }
}


void State::RegisterDomainSet(const Key& name,
        const Teuchos::RCP<AmanziMesh::DomainSet> set) {
  domain_sets_[name] = set;
}

bool State::HasDomainSet(const Key& name) const {
  return (bool) domain_sets_.count(name);
}

Teuchos::RCP<const AmanziMesh::DomainSet>
State::GetDomainSet(const Key& name) const
{
  if (!domain_sets_.count(name)) {
    Errors::Message msg;
    msg << "DomainSet \"" << name << "\" does not exist in State.";
    Exceptions::amanzi_throw(msg);
  }
  return domain_sets_.at(name);
}


// -----------------------------------------------------------------------------
// State handles data evaluation.
// -----------------------------------------------------------------------------
Evaluator& State::RequireEvaluator(const Key& key, const Key& tag)
{
  // does it already exist?
  if (HasEvaluator(key, tag)) {
    return GetEvaluator(key, tag);
  }

  // See if the key is provided by another existing evaluator.
  for (auto& e : evaluators_) {
    if (Keys::hasKey(e.second, tag) &&
        e.second.at(tag)->ProvidesKey(key, tag)) {
      auto& evaluator = e.second.at(tag);
      SetEvaluator(key, evaluator);
      return *evaluator;
    }
  }

  // Create the evaluator from State's plist
  // -- Get the Field Evaluator plist
  Teuchos::ParameterList& fm_plist = state_plist_.sublist("evaluators");

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
      Teuchos::Array<std::string> modelnames = sublist.get<Teuchos::Array<std::string>>("models parameters");
      for (auto modelname = modelnames.begin(); modelname != modelnames.end(); ++modelname) {
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
    return *evaluator;
  }

  // cannot find the evaluator, error
  Errors::Message message;
  message << "Model for \"" << key << "\" cannot be created in State.";
  Exceptions::amanzi_throw(message);
}


bool State::HasEvaluator(const Key& key, const Key& tag)
{
  if (Keys::hasKey(evaluators_, key)) {
    return Keys::hasKey(evaluators_.at(key), tag);
  } else {
    return false;
  }
}


Teuchos::RCP<const Functions::MeshPartition> State::GetMeshPartition(Key key)
{
  Teuchos::RCP<const Functions::MeshPartition> mp = GetMeshPartition_(key);
  if (mp == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Mesh partition " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return mp;
}


Teuchos::RCP<const Functions::MeshPartition> State::GetMeshPartition_(Key key)
{
  MeshPartitionMap::iterator lb = mesh_partitions_.lower_bound(key);
  if (lb != mesh_partitions_.end() && !(mesh_partitions_.key_comp()(key, lb->first))) {
    return lb->second;
  } else {
    if (state_plist_.isParameter("mesh partitions")) {
      Teuchos::ParameterList& part_superlist = state_plist_.sublist("mesh partitions");
      if (part_superlist.isParameter(key)) {
        Teuchos::ParameterList& part_list = part_superlist.sublist(key);
        std::string mesh_name = part_list.get<std::string>("mesh name", "domain");
        Teuchos::Array<std::string> region_list =
            part_list.get<Teuchos::Array<std::string> >("region list");
        Teuchos::RCP<Functions::MeshPartition> mp = Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::CELL,
                region_list.toVector()));
        mp->Initialize(GetMesh(mesh_name), -1);
        mp->Verify();
        mesh_partitions_[key] = mp;
        return mp;
      }
    }
    return Teuchos::null;
  }
}


void State::WriteDependencyGraph() const
{
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    std::ofstream os("dependency_graph.txt", std::ios::out);
    for (auto& e : evaluators_) {
      for (auto& r : e.second) os << *r.second;
    }
    os.close();
  }
}


// -----------------------------------------------------------------------------
// State handles model parameters.
// -----------------------------------------------------------------------------
Teuchos::ParameterList
State::GetModelParameters(std::string modelname) {
  AMANZI_ASSERT(state_plist_.isSublist("model parameters"));
  Teuchos::ParameterList model_plist = state_plist_.sublist("model parameters");
  AMANZI_ASSERT(model_plist.isSublist(modelname));
  return model_plist.sublist(modelname);
}


// -----------------------------------------------------------------------------
// Initialization, etc.
// -----------------------------------------------------------------------------
// Initialize data, allowing values to be specified here or in the owning PK.
// All independent variables must be initialized here.
void State::Setup()
{
  Require<double>("time", StateTags::DEFAULT, "time");
  Require<double>("time", "next", "time");
  Require<int>("cycle", StateTags::DEFAULT, "cycle");
  Require<int>("cycle", "next", "cycle");

  Teuchos::OSTab tab = vo_->getOSTab();

  // Ensure compatibility of all the evaluators -- each evaluator's dependencies
  // must provide what is required of that evaluator.
  for (auto& e : evaluators_) {
    for (auto& r : e.second) {
      if (!r.second->ProvidesKey(e.first, r.first)) {
        Errors::Message msg;
        msg << "Evaluator \"" << e.first << "\" with tag \"" << r.first
            << "\" does not provide its own key.";
        Exceptions::amanzi_throw(msg);
      }

      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        Teuchos::OSTab tab1 = vo_->getOSTab();
        *vo_->os() << "checking compatibility: " << r.first;
      }
      r.second->EnsureCompatibility(*this);
    }
  }

  // -- Create the data for all fields.
  for (auto& r : data_) {
    r.second->CreateData();

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Record \"" << r.first << "\": ";

      auto& field = GetRecord(r.first);
      if (field.ValidType<CompositeVector>()) {
        const auto& cv = Get<CompositeVector>(r.first);
        for (auto comp = cv.begin(); comp != cv.end(); ++comp) *vo_->os() << *comp << " ";
      } else {
        *vo_->os() << " not of type CV";
      }
      *vo_->os() << "\n";
    }
  }

  for (auto& deriv : derivs_) {
    deriv.second->CreateData();
  }
}


void State::Initialize()
{
  // Set metadata
  GetW<int>("cycle", "", "cycle") = cycle_;
  data_.at("cycle")->initializeTags();

  GetW<double>("time", "", "time") = time_;
  data_.at("time")->initializeTags();

  // Initialize data from initial conditions.
  InitializeFields();

  // Initialize evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  CheckAllFieldsInitialized();

  // Write dependency graph.
  WriteDependencyGraph();

  // Reset io_vis flags using blacklist and whitelist
  InitializeIOFlags();
}


void State::Initialize(Teuchos::RCP<State> S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab1 = vo_->getOSTab();
    *vo_->os() << "copying fields to new state.." << std::endl;
  }

  for (auto& e : data_) {
    if (S->HasData(e.first)) {
      auto& field = GetRecordW(e.first, e.first);
      const auto& copy = S->GetRecord(e.first);

      /*
      if (field.type() != copy.type()) {
        std::stringstream messagestream;
        messagestream << "States has fields with the same name but different types\n";
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
      }
      */
      if (copy.initialized()) {
        if (field.ValidType<double>()) {
          field.Set<double>(field.owner(), copy.Get<double>());
        } else if (field.ValidType<std::vector<double> >()) {
          field.Set<std::vector<double>>(field.owner(), copy.Get<std::vector<double>>());
        } else if (field.ValidType<CompositeVector>()) {
          field.Set<CompositeVector>(field.owner(), copy.Get<CompositeVector>());
        } else {
          Errors::Message message("Copy field with unknown type\n");
          Exceptions::amanzi_throw(message);
        }
        field.set_initialized();      
      }
    }   
  }

  // Initialize any other fields from state plist.
  InitializeFields();

  // Initialize other field evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  // CheckAllFieldsInitialized();

  // Write dependency graph.
  WriteDependencyGraph();

  // Reset io_vis flags using blacklist and whitelist
  InitializeIOFlags();
}


void State::InitializeEvaluators()
{
  for (auto& e : evaluators_) {
    auto tag = e.second.begin();

    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "initializing eval: \"" << e.first << "\" with 1st tag \"" << tag->first << "\"" << std::endl;
    }

    tag->second->Update(*this, "state");
    GetRecordSetW(e.first).initializeTags();
  }
}


void State::InitializeFieldCopies()
{
  VerboseObject vo("State", state_plist_); 
  if (vo.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "initializing field copies..." << std::endl;
  }

  for (auto& e : data_) {
    const auto& copy = GetRecord(e.first);
    if (copy.ValidType<CompositeVector>()) {
      for (auto& r : *e.second) {
        Set<CompositeVector>(e.first, r.first, copy.owner(), copy.Get<CompositeVector>());
      }
    } 
  }
}


void State::InitializeFields()
{
  bool pre_initialization = false;
  VerboseObject vo("State", state_plist_); 

  if (state_plist_.isParameter("initialization filename")) {
    pre_initialization = true;
    std::string filename = state_plist_.get<std::string>("initialization filename");
    Amanzi::Checkpoint file_input(filename, GetMesh()->get_comm());
    // file_input.open_h5file();

    for (auto it = data_.begin(); it != data_.end(); ++it) {
      auto owner = GetRecord(it->first).owner();
      auto& r = GetRecordW(it->first, owner); 
      if (r.ValidType<CompositeVector>()) {
        r.ReadCheckpoint(file_input);

        if (HasEvaluator(it->first)) {
          Evaluator& fm = GetEvaluator(it->first);
          auto tmp = dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>* >(&fm);
          if (tmp != nullptr) {
            tmp->SetChanged();
          }
        }          
        it->second->initializeTags();
      }
    }
    // file_input.close_h5file();
  }
   
  // Initialize through initial condition
  if (vo.os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "initializing data through initial conditions..." << std::endl;
  }

  if (state_plist_.isSublist("initial conditions")) {
    for (auto& e : data_) {
      bool flag(false);
      if (pre_initialization || !e.second->isInitialized()) {
        if (state_plist_.sublist("initial conditions").isSublist(e.first)) {
          flag = true;
          Teuchos::ParameterList sublist = state_plist_.sublist("initial conditions").sublist(e.first);
          e.second->Initialize(sublist);
        } else {
          // check for domain set
          KeyTriple split;
          bool is_ds = Keys::splitDomainSet(e.first, split);
          Key ds_name = std::get<0>(split);
          if (is_ds) {
            Key lifted_key = Keys::getKey(ds_name, "*", std::get<2>(split));
            if (state_plist_.sublist("initial conditions").isSublist(lifted_key)) {
              flag = true;
              Teuchos::ParameterList sublist = state_plist_.sublist("initial conditions").sublist(lifted_key);
              sublist.set("evaluator name", e.first);
              e.second->Initialize(sublist);
            }
          }
        }
      }

      if (vo.os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo.getOSTab();
        *vo_->os() << "initializing \"" << e.first << "\"  [" << flag << "]" << std::endl;
      }
    }
  }
}


// Make sure all fields have gotten their IC, either from State or the owning PK.
bool State::CheckAllFieldsInitialized()
{
  for (auto& e : data_) {
    if (!e.second->isInitialized()) {
      // field was not initialized
      std::stringstream ss;
      ss << "Variable \"" << e.first << "\" was not initialized\n";
      Errors::Message msg(ss.str());
      Exceptions::amanzi_throw(msg);
      return false;
    }
  }
  return true;
}


// Utility for setting vis flags
void State::InitializeIOFlags()
{
  Teuchos::Array<std::string> empty;

  // removing fields from vis dump
  std::vector<std::string> blacklist = 
      state_plist_.get<Teuchos::Array<std::string> >("blacklist", empty).toVector();

  for (auto it = data_begin(); it != data_end(); ++it) {
    bool io_block(false);
    for (int m = 0; m < blacklist.size(); ++m) {
      std::regex pattern(blacklist[m]);
      io_block |= std::regex_match(it->first, pattern);
    }
    for (auto& r : *it->second) r.second->set_io_vis(!io_block);
  }

  // adding fields to vis dump
  std::vector<std::string> whitelist =
      state_plist_.get<Teuchos::Array<std::string> >("whitelist", empty).toVector();

  for (auto it = data_begin(); it != data_end(); ++it) {
    bool io_allow(false);
    for (int m = 0; m < whitelist.size(); ++m) {
      std::regex pattern(whitelist[m]);
      io_allow |= std::regex_match(it->first, pattern);
    }
    if (io_allow) { 
      for (auto& r : *it->second) r.second->set_io_vis(true);
    }
  }
}


// Non-member function for vis.
void WriteVis(Visualization& vis, const State& S)
{
  if (!vis.is_disabled()) {
    // Create the new time step
    vis.CreateTimestep(S.time(), S.cycle(), vis.get_tag());

    // Write all fields to the visualization file, the fields know if they
    // need to be written.
    for (auto r = S.data_begin(); r != S.data_end(); ++r) {
      r->second->WriteVis(vis);
    }
    vis.WriteRegions();
    vis.WritePartition();
    vis.FinalizeTimestep();
  }
}


// Non-member function for checkpointing.
void WriteCheckpoint(Checkpoint& chkp, const Comm_ptr_type& comm,
                     const State& S, bool final)
{
  if (!chkp.is_disabled()) {
    // chkp.SetFinal(final);
    chkp.CreateFile(S.cycle());
    for (auto r = S.data_begin(); r != S.data_end(); ++r) {
      r->second->WriteCheckpoint(chkp);
    }
    chkp.Write("mpi_num_procs", comm->NumProc());
    chkp.Finalize();
  }
}


// Non-member function for checkpointing.
void ReadCheckpoint(const Comm_ptr_type& comm, State& S, 
                    const std::string& filename)
{
  Checkpoint chkp(filename, comm);

  // Load the number of processes and ensure they are the same.
  int num_procs(-1);
  chkp.Read("mpi_num_procs", num_procs);
  if (comm->NumProc() != num_procs) {
    std::stringstream ss;
    ss << "Requested checkpoint file was created on " << num_procs
       << " processes, making it incompatible with this run on "
       << comm->NumProc() << " processes.";
    Errors::Message message(ss.str());
    throw(message);
  }

  // load the data
  for (auto data = S.data_begin(); data != S.data_end(); ++data) {
    data->second->ReadCheckpoint(chkp);
  }

  chkp.Finalize();
}


// Non-member function for checkpointing.
double ReadCheckpointInitialTime(const Comm_ptr_type& comm,
                                 std::string filename)
{
  if (!Keys::ends_with(filename, ".h5")) {
    // new style checkpoint
    boost::filesystem::path filepath = boost::filesystem::path(filename) / "domain.h5";
    filename = filepath.string();
  }
  HDF5_MPI checkpoint(comm, filename);

  // load the attributes
  double time(0.);
  checkpoint.open_h5file();
  checkpoint.readAttrReal(time, "time");
  checkpoint.close_h5file();
  return time;
}


// Non-member function for checkpointing position.
int ReadCheckpointPosition(const Comm_ptr_type& comm, std::string filename)
{
  if (!Keys::ends_with(filename, ".h5")) {
    // new style checkpoint
    boost::filesystem::path filepath = boost::filesystem::path(filename) / "domain.h5";
    filename = filepath.string();
  }
  HDF5_MPI checkpoint(comm, filename);

  // load the attributes
  int pos = 0;
  checkpoint.open_h5file();
  checkpoint.readAttrInt(pos, "position");
  checkpoint.close_h5file();
  return pos;
}


// Non-member function for checkpointing observations.
void ReadCheckpointObservations(const Comm_ptr_type& comm,
                                std::string filename,
                                Amanzi::ObservationData& obs_data)
{
  if (!Keys::ends_with(filename, ".h5")) {
    // new style checkpoint
    boost::filesystem::path filepath = boost::filesystem::path(filename) / "domain.h5";
    filename = filepath.string();
  }

  HDF5_MPI checkpoint(comm, filename);
  checkpoint.open_h5file();

  // read observations
  int nlabels, ndata(0), ndata_glb(0);
  int* nobs;
  char** tmp_labels;
  double* tmp_data(NULL);

  checkpoint.readDataString(&tmp_labels, &nlabels, "obs_names");
  if (nlabels > 0) { 
    checkpoint.readAttrInt(&nobs, &nlabels, "obs_numbers");
  }
  for (int i = 0; i < nlabels; ++i) ndata_glb += 2 * nobs[i];
  ndata = (comm->MyPID() == 0) ? ndata_glb : 0;
  checkpoint.readDatasetReal(&tmp_data, ndata, "obs_values");

  checkpoint.close_h5file();

  // populated observations on root
  if (comm->MyPID() == 0) {
    int m(0);
    Amanzi::ObservationData::DataQuadruple data_quad;

    for (int i = 0; i < nlabels; ++i) {
      std::vector<ObservationData::DataQuadruple>& od = obs_data[tmp_labels[i]];
      for (int k = 0; k < nobs[i]; ++k) {
        data_quad.time = tmp_data[m++];
        data_quad.value = tmp_data[m++];
        data_quad.is_valid = true;
        od.push_back(data_quad);
      }
    }
  }

  // clean memory
  for (int i = 0; i < nlabels; i++) free(tmp_labels[i]);
  if (nlabels > 0) {
    free(tmp_labels);
    free(nobs);
    if (tmp_data != NULL) free(tmp_data); 
  }
}


// Non-member function for deforming the mesh after reading a checkpoint file
// that contains the vertex coordinate field (this is written by deformation pks)
//
// FIX ME: Refactor this to make the name more general.  Should align with a
// mesh name prefix or something, and the coordinates should be written by
// state in WriteCheckpoint if mesh IsDeformableMesh() --ETC
void DeformCheckpointMesh(State& S, Key domain)
{
  if (S.HasData("vertex coordinate")) { // only deform mesh if vertex coordinate
                                        // field exists
    AmanziMesh::Mesh* write_access_mesh = const_cast<AmanziMesh::Mesh*>(&*S.GetMesh());

    // get vertex coordinates state
    const CompositeVector& vc = S.Get<CompositeVector>("vertex coordinate");
    vc.ScatterMasterToGhosted("node");
    const Epetra_MultiVector& vc_n = *vc.ViewComponent("node", true);

    int dim = write_access_mesh->space_dimension();
    Amanzi::AmanziMesh::Entity_ID_List nodeids;
    Amanzi::AmanziGeometry::Point new_coords(dim);
    AmanziGeometry::Point_List new_pos, final_pos;

    int nV = vc_n.MyLength();
    for (int n = 0; n != nV; ++n) {
      for (int k = 0; k != dim; ++k)
        new_coords[k] = vc_n[k][n];

      // push back for deform method
      nodeids.push_back(n);
      new_pos.push_back(new_coords);
    }

    // deform the mesh
    if (Keys::starts_with(domain, "column"))
      write_access_mesh->deform(nodeids, new_pos, false, &final_pos);
    else
      write_access_mesh->deform(nodeids, new_pos, true, &final_pos);
  }
}


// Non-member function for statistics
void WriteStateStatistics(const State& S, const VerboseObject& vo)
{
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "\nField                                    Min/Max/Avg" << std::endl;

    for (auto it = S.data_begin(); it != S.data_end(); ++it) {
      std::string name(it->first);

      if (name.size() > 33) replace_all(name, "temperature", "temp");
      if (name.size() > 33) replace_all(name, "internal_energy", "ie");
      if (name.size() > 33) replace_all(name, "molar", "mol");

      if (S.GetRecord(name).ValidType<CompositeVector>()) {
        std::map<std::string, double> vmin, vmax, vavg;
        S.Get<CompositeVector>(name).MinValue(vmin);
        S.Get<CompositeVector>(name).MaxValue(vmax);
        S.Get<CompositeVector>(name).MeanValue(vavg);

        for (auto c_it = vmin.begin(); c_it != vmin.end(); ++c_it) {
          std::string namedot(name), name_comp(c_it->first);
          if (vmin.size() != 1) namedot.append("." + name_comp);
          namedot.resize(40, '.');
          *vo.os() << namedot << " " << c_it->second << " / " 
                    << vmax[name_comp] << " / " << vavg[name_comp] << std::endl;
        }
      }
      else if (S.GetRecord(name).ValidType<double>()) {
        double vmin = S.Get<double>(name);
        name.resize(40, '.');
        *vo.os() << name << " " << vmin << std::endl;
      }
      else if (S.GetRecord(name).ValidType<AmanziGeometry::Point>()) {
        const auto& p = S.Get<AmanziGeometry::Point>(name);
        name.resize(40, '.');
        *vo.os() << name;
        for (int i = 0; i < p.dim(); ++i) *vo.os() << " " << p[i];
        *vo.os() << std::endl;
      }
    }
  }
}


Evaluator& State::GetEvaluator(const Key& key, const Key& tag) {
  return *GetEvaluatorPtr(key, tag);
}


const Evaluator& State::GetEvaluator(const Key& key, const Key& tag) const
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
}


Teuchos::RCP<Evaluator> State::GetEvaluatorPtr(const Key& key, const Key& tag)
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
}


void State::SetEvaluator(const Key& key, const Key& tag,
                         const Teuchos::RCP<Evaluator>& evaluator) {
  evaluators_[key][tag] = evaluator;
}


Teuchos::ParameterList&
State::GetEvaluatorList(const Key& key)
{
  if (FEList().isParameter(key)) {
    return FEList().sublist(key);
  } else {
    // check for domain set
    KeyTriple split;
    bool is_ds = Keys::splitDomainSet(key, split);
    if (is_ds) {
      Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
      if (FEList().isParameter(lifted_key)) {
        return FEList().sublist(lifted_key);
      }
    }
  }

  return FEList().sublist(key);
}

}  // namespace Amanzi
