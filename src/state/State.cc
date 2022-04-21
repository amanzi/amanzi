/*
  State

  Copyright 2010-202x held jointly by LANS/LANL, LBNL, and PNNL.
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

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------
State::State() {
  vo_ = Teuchos::rcp(new VerboseObject("State", "none"));
};

State::State(Teuchos::ParameterList& state_plist) :
    state_plist_(state_plist) {
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


void State::AliasMesh(const Key& target, const Key& alias) {
  bool deformable = IsDeformableMesh(target);
  Teuchos::RCP<AmanziMesh::Mesh> mesh = GetMesh_(target);
  RegisterMesh(alias, mesh, deformable);
  mesh_aliases_[alias] = target;

  if (GetMesh_(target+"_3d") != Teuchos::null) {
    RegisterMesh(alias+"_3d", GetMesh_(target+"_3d"), deformable);
    mesh_aliases_[alias+"_3d"] = target+"_3d";
  }
};


bool State::IsAliasedMesh(const Key& key) const {
  return Keys::hasKey(mesh_aliases_, key);
}


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
State::GetDomainSet(const Key& name) const {
  if (!domain_sets_.count(name)) {
    Errors::Message msg;
    msg << "DomainSet \"" << name << "\" does not exist in State.";
    Exceptions::amanzi_throw(msg);
  }
  return domain_sets_.at(name);
}

// -----------------------------------------------------------------------------
// State handles data
// -----------------------------------------------------------------------------
//
// RecordSets
//
bool State::HasRecordSet(const Key& fieldname) const {
  return Keys::hasKey(data_, fieldname);
}

const RecordSet& State::GetRecordSet(const Key& fieldname) const {
  if (!HasRecordSet(fieldname)) {
    Errors::Message msg;
    msg << "State does not have a RecordSet named \"" << fieldname << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return *data_.at(fieldname);
}

RecordSet& State::GetRecordSetW(const Key& fieldname) {
  if (!HasRecordSet(fieldname)) {
    Errors::Message msg;
    msg << "State does not have a RecordSet named \"" << fieldname << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return *data_.at(fieldname);
}


//
// Records
//
bool State::HasRecord(const Key& fieldname, const Tag& tag) const {
  return HasRecordSet(fieldname) && GetRecordSet(fieldname).HasRecord(tag);
}

const Record& State::GetRecord(const Key& fieldname, const Tag& tag) const {
  return GetRecordSet(fieldname).GetRecord(tag);
}

Record& State::GetRecordW(const Key& fieldname, const Key& owner) {
  auto& r = GetRecordSetW(fieldname).GetRecord(Tags::DEFAULT);
  r.AssertOwnerOrDie(owner);
  return r;
}
Record& State::GetRecordW(const Key& fieldname, const Tag& tag, const Key& owner) {
  auto& r = GetRecordSetW(fieldname).GetRecord(tag);
  r.AssertOwnerOrDie(owner);
  return r;
}

//
// Derivative RecordSets
//
bool State::HasDerivativeSet(const Key& fieldname, const Tag& tag) const {
  return Keys::hasKey(derivs_, Keys::getKey(fieldname, tag));
}

const RecordSet& State::GetDerivativeSet(const Key& fieldname, const Tag& tag) const {
  if (!HasDerivativeSet(fieldname, tag)) {
    Errors::Message msg;
    msg << "State does not have a Derivative RecordSet for field \"" << fieldname << "\" at tag \""
        << tag.get() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return *derivs_.at(Keys::getKey(fieldname, tag));
}

RecordSet& State::GetDerivativeSetW(const Key& fieldname, const Tag& tag) {
  if (!HasDerivativeSet(fieldname, tag)) {
    Errors::Message msg;
    msg << "State does not have a Derivative RecordSet for field \"" << fieldname << "\" at tag \""
        << tag.get() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return *derivs_.at(Keys::getKey(fieldname, tag));
}

//
// Derivative Records -- note, there is currently no interface for directly
// getting the Derivative Record.  We expect this is not necessary.
//
bool State::HasDerivative(const Key& key, const Tag& tag,
                          const Key& wrt_key, const Tag& wrt_tag) const {
  auto keytag = Keys::getKey(key, tag);
  if (Keys::hasKey(derivs_, keytag)) {
    Tag der_tag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    return derivs_.at(keytag)->HasRecord(der_tag);
  }
  return false;
}
bool State::HasDerivative(const Key& key, const Key& wrt_key) const {
  return HasDerivative(key, Tags::DEFAULT, wrt_key, Tags::DEFAULT);
}


// -----------------------------------------------------------------------------
// State handles data evaluation.
// -----------------------------------------------------------------------------
Evaluator& State::RequireEvaluator(const Key& key, const Tag& tag)
{
  CheckIsDebugEval_(key, tag);

  // does it already exist?
  if (HasEvaluator(key, tag)) {
    return GetEvaluator(key, tag);
  }

  // See if the key is provided by another existing evaluator.
  for (auto& e : evaluators_) {
    if (Keys::hasKey(e.second, tag) &&
        e.second.at(tag)->ProvidesKey(key, tag)) {
      auto& evaluator = e.second.at(tag);
      SetEvaluator(key, tag, evaluator);
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
    sublist.set("tag", tag.get());
    auto evaluator = evaluator_factory.createEvaluator(sublist);
    SetEvaluator(key, tag, evaluator);
    return *evaluator;
  }

  // is it cell volume?
  if (Keys::getVarName(key) == "cell_volume") {
    Teuchos::ParameterList& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "cell volume");
    return RequireEvaluator(key, tag);
  }

  // cannot find the evaluator, error
  Errors::Message message;
  message << "Model for \"" << key << "@" << tag.get() << "\" cannot be created in State.";
  Exceptions::amanzi_throw(message);
  return *Evaluator_Factory().createEvaluator(fm_plist); // silences warning
}


bool State::HasEvaluator(const Key& key, const Tag& tag)
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
  // FIXME -- this is not what it used to be.  This simply writes data
  // struture, and is not the dependency graph information at all.  Rename
  // this, then recover the old WriteDependencyGraph method, which wrote a list
  // of all evaluators and their dependnecies that could be read in networkx
  // for plotting the dag. --ETC

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "------------------------------------------" << std::endl
               << "Dependency & Structure list for evaluators" << std::endl
               << "------------------------------------------" << std::endl;
    for (auto& e : evaluators_) {
      for (auto& r : e.second) *vo_->os() << *r.second;
      if (GetRecord(e.first, e.second.begin()->first).ValidType<CompositeVector>()) {
        Teuchos::OSTab tab1 = vo_->getOSTab();
        Teuchos::OSTab tab2 = vo_->getOSTab();
        data_.at(e.first)->GetFactory<CompositeVector, CompositeVectorSpace>().Print(*vo_->os());
        *vo_->os() << std::endl;
      }
    }

    *vo_->os() << "------------------------------" << std::endl
               << "Structure list for derivatives" << std::endl
               << "------------------------------" << std::endl;
    for (auto& e : derivs_) {
      *vo_->os() << "D" << e.first << "/D{ ";
      for (const auto& wrt : *e.second) *vo_->os() << wrt.first.get() << " ";
      *vo_->os() << "}" << std::endl;
      auto wrt_tag = e.second->begin();
      if (e.second->GetRecord(wrt_tag->first).ValidType<CompositeVector>()) {
        Teuchos::OSTab tab1 = vo_->getOSTab();
        Teuchos::OSTab tab2 = vo_->getOSTab();
        e.second->GetFactory<CompositeVector, CompositeVectorSpace>().Print(*vo_->os());
        *vo_->os() << std::endl;
      }
    }
  }
}


// -----------------------------------------------------------------------------
// State handles model parameters.
// -----------------------------------------------------------------------------
Teuchos::ParameterList
State::GetModelParameters(std::string modelname)
{
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
  require_time(Tags::DEFAULT);
  GetRecordSetW("time").initializeTags();

  require_cycle(Tags::DEFAULT);
  GetRecordSetW("cycle").initializeTags();

  Require<double>("dt", Tags::DEFAULT, "dt");
  GetRecordW("dt", Tags::DEFAULT, "dt").set_initialized();

  Require<int>("position", Tags::DEFAULT, "position");
  GetRecordSetW("position").initializeTags();

  Teuchos::OSTab tab = vo_->getOSTab();

  // Ensure consistency of the dependency graph.  This takes two steps -- the
  // first pass ensures that all evaluators are present and instantiated,
  // completed the dag.  The second pass sets data requirements, and makes sure
  // data requirements are consistent.
  //
  // Note that the first pass may modify the graph, but since it is a DAG, and
  // this is called recursively, we can just call it on the nodes that appear
  // initially.
  { // scope for copy
    EvaluatorMap evaluators_copy(evaluators_);
    for (auto& e : evaluators_copy) {
      for (auto& r : e.second) {
        // if (!r.second->ProvidesKey(e.first, r.first)) {
        //   Errors::Message msg;
        //   msg << "Evaluator \"" << e.first << "\" with tag \"" << r.first.get()
        //       << "\" does not provide its own key.";
        //   Exceptions::amanzi_throw(msg);
        // }

        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          Teuchos::OSTab tab1 = vo_->getOSTab();
          *vo_->os() << "ensure evaluators: \"" << e.first << "\" @ \""
                     << r.first << "\"" << std::endl;
        }
        r.second->EnsureEvaluators(*this);
      }
    }
  }

  // Second pass calls EnsureCompatibility, which checks data consistency.
  // This pass does not modify the graph.
  for (auto& e : evaluators_) {
    for (auto& r : e.second) {
      r.second->EnsureCompatibility(*this);
    }
  }

  // Create the data for all fields.
  for (auto& r : data_) {
    r.second->CreateData();

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "RecordSet \"" << r.first << "\", tags: ";

      for(auto& e : *r.second) *vo_->os() << "\"" << e.first.get() << "\" ";
      if (r.second->ValidType<CompositeVector, CompositeVectorSpace>()) {
        const auto& cvs = r.second->GetFactory<CompositeVector, CompositeVectorSpace>();
        *vo_->os() << "comps: ";
        for (const auto& comp : cvs) *vo_->os() << comp << " ";
      } else {
        *vo_->os() << "not CV";
      }
      *vo_->os() << "\n";
    }
  }

  // Create the data for all derivatives
  for (auto& deriv : derivs_) {
    deriv.second->CreateData();
  }

  // -- Write DAG to disk for visualization
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    WriteDependencyGraph();
    *vo_->os() << "Setup is complete.\n\n";
  }
}


void State::Initialize()
{
  // Set metadata
  GetW<int>("cycle", Tags::DEFAULT, "cycle") = 0;
  GetRecordSetW("cycle").initializeTags();

  GetW<double>("time", Tags::DEFAULT, "time") = 0.0;
  GetRecordSetW("time").initializeTags();

  GetW<double>("dt", Tags::DEFAULT, "dt") = 0.;
  GetRecordSetW("dt").initializeTags();

  // Initialize data from initial conditions.
  InitializeFields();

  // Initialize evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  CheckAllFieldsInitialized();

  // Reset io_vis flags using blacklist and whitelist
  InitializeIOFlags();
}


//
// NOTE: this method assumes all Records have a DEFAULT tag
//
void State::Initialize(const State& other)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab1 = vo_->getOSTab();
    *vo_->os() << "copying fields to new state.." << std::endl;
  }

  for (auto& e : data_) {
    if (other.HasRecord(e.first, Tags::DEFAULT)) {
      auto owner = GetRecord(e.first, Tags::DEFAULT).owner();
      auto& field = GetRecordW(e.first, Tags::DEFAULT, owner);
      const auto& copy = other.GetRecord(e.first, Tags::DEFAULT);

      if (copy.initialized()) {
        field.Assign(copy);
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

  // Reset io_vis flags using blacklist and whitelist
  InitializeIOFlags();
}


void State::InitializeEvaluators()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  for (const auto& e : evaluators_) {
    for (const auto& tag : e.second) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "initializing eval: \"" << e.first << "\" @ \"" << tag.first << "\"" << std::endl;
      }

      tag.second->Update(*this, "state");
    }
    GetRecordSetW(e.first).initializeTags();
  }
}


void State::InitializeFieldCopies(const Tag& reference_tag)
{
  VerboseObject vo("State", state_plist_);
  if (vo.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "initializing field copies..." << std::endl;
  }

  for (auto& e : data_) {
    if (HasRecord(e.first, reference_tag)) {
      const auto& copy = GetRecord(e.first, reference_tag);
      if (copy.initialized() && copy.ValidType<CompositeVector>()) {
        for (auto& r : *e.second) {
          if (!r.second->initialized()) {
            r.second->Assign(copy);
            r.second->set_initialized();
          }
        }
      }
    }
  }
}


void State::InitializeFields(const Tag& tag)
{
  bool pre_initialization = false;
  VerboseObject vo("State", state_plist_);

  // this directly overlaps with "initial conditions" --> "varname" -->
  // "restart file" capability, but is done globally instead of by variable.
  // Can the be refactored to use that instead?
  if (state_plist_.isParameter("initialization filename")) {
    pre_initialization = true;
    std::string filename = state_plist_.get<std::string>("initialization filename");
    Amanzi::Checkpoint file_input(filename, GetMesh()->get_comm());
    // file_input.open_h5file();

    for (auto it = data_.begin(); it != data_.end(); ++it) {
      auto owner = GetRecord(it->first, tag).owner();
      auto& r = GetRecordW(it->first, tag, owner);
      if (r.ValidType<CompositeVector>()) {
        r.ReadCheckpoint(file_input);

        // this is pretty hacky -- why are these ICs not in the PK's list?  And
        // if they aren't owned by a PK, they should be independent variables
        // not primary variables.
        if (HasEvaluator(it->first, tag)) {
          Evaluator& fm = GetEvaluator(it->first, tag);
          auto tmp = dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>* >(&fm);
          if (tmp != nullptr) {
            tmp->SetChanged();
          }
        }
        it->second->initializeTags();
      }

      if (vo.os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo.getOSTab();
        *vo_->os() << "initializing from HDF5 file: \"" << it->first << std::endl;
      }
    }
    // file_input.close_h5file();
  }

  // Initialize through initial condition
  if (vo.os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "initializing data through initial conditions..." << std::endl;
  }

  Tag failed;
  if (state_plist_.isSublist("initial conditions")) {
    for (auto& e : data_) {
      std::string flag("... skipped");
      if (pre_initialization || !e.second->isInitialized(failed)) {
        if (state_plist_.sublist("initial conditions").isSublist(e.first)) {
          flag = "[ok]";
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
              flag = "[ok]";
              Teuchos::ParameterList sublist = state_plist_.sublist("initial conditions").sublist(lifted_key);
              sublist.set("evaluator name", e.first);
              e.second->Initialize(sublist);
            }
          }
        }
      }

      if (vo.os_OK(Teuchos::VERB_HIGH)) {
        Teuchos::OSTab tab = vo.getOSTab();
        *vo_->os() << "initializing \"" << e.first << "\" " << flag << std::endl;
      }
    }
  }
}


// Make sure all fields have gotten their IC, either from State or the owning PK.
bool State::CheckAllFieldsInitialized()
{
  Tag failed;
  for (auto& e : data_) {
    if (!e.second->isInitialized(failed)) {
      std::stringstream ss;
      ss << "Variable \"" << e.first << "\" with tag \"" << failed.get() << "\" was not initialized\n";
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


Evaluator& State::GetEvaluator(const Key& key, const Tag& tag) {
  return *GetEvaluatorPtr(key, tag);
}


const Evaluator& State::GetEvaluator(const Key& key, const Tag& tag) const
{
  try {
    return *evaluators_.at(key).at(tag);
  } catch (std::out_of_range) {
    std::stringstream ss;
    ss << "Evaluator for field \"" << key << "\" at tag \"" << tag
       << "\" does not exist in the state.";
    Errors::Message message(ss.str());
    throw(message);
  }
}


Teuchos::RCP<Evaluator> State::GetEvaluatorPtr(const Key& key, const Tag& tag)
{
  try {
    return evaluators_.at(key).at(tag);
  } catch (std::out_of_range) {
    std::stringstream ss;
    ss << "Evaluator for field \"" << key << "\" at tag \"" << tag
       << "\" does not exist in the state.";
    Errors::Message message(ss.str());
    throw(message);
  }
}


// Key defines field, Tag defines particular copy of field
void State::SetEvaluator(const Key& key, const Tag& tag,
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

  // return an empty new list
  return FEList().sublist(key);
}


bool
State::HasEvaluatorList(const Key& key)
{
  if (FEList().isSublist(key)) return true;
  // check for domain set
  KeyTriple split;
  bool is_ds = Keys::splitDomainSet(key, split);
  if (is_ds) {
    Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
    if (FEList().isSublist(lifted_key)) return true;
  }
  return false;
}


void State::CheckIsDebugEval_(const Key& key, const Tag& tag)
{
  // check for debugging.  This provides a line for setting breakpoints for
  // debugging PK and Evaluator dependencies.
#ifdef ENABLE_DBC
  Teuchos::Array<Key> debug_evals = state_plist_.sublist("debug").get<Teuchos::Array<std::string>>("evaluators",
          Teuchos::Array<Key>());
  if (std::find(debug_evals.begin(), debug_evals.end(), key) != debug_evals.end()) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "State: Evaluator for debug field \"" << key << "@" << tag << "\" was required." << std::endl;
    }
  }
#endif
}

void State::CheckIsDebugData_(const Key& key, const Tag& tag)
{
  // check for debugging.  This provides a line for setting breakpoints for
  // debugging PK and Evaluator dependencies.
#ifdef ENABLE_DBC
  Teuchos::Array<Key> debug_evals = state_plist_.sublist("debug").get<Teuchos::Array<std::string>>("data",
          Teuchos::Array<Key>());
  if (std::find(debug_evals.begin(), debug_evals.end(), key) != debug_evals.end()) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "State: data for debug field \"" << key << "@" << tag << "\" was required." << std::endl;
    }
  }
#endif
}

}  // namespace Amanzi
