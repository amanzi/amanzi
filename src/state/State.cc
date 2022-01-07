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
Evaluator& State::RequireEvaluator(const Key& key, const Tag& tag)
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
    sublist.set("tag", tag.get());
    auto evaluator = evaluator_factory.createEvaluator(sublist);
    SetEvaluator(key, tag, evaluator);
    return *evaluator;
  }

  // cannot find the evaluator, error
  Errors::Message message;
  message << "Model for \"" << key << "\" cannot be created in State.";
  Exceptions::amanzi_throw(message);
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
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    for (auto& e : evaluators_) {
      for (auto& r : e.second) *vo_->os() << *r.second;
    }
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
  Require<double>("time", Tags::DEFAULT, "time");
  Require<double>("time", Tags::NEXT, "time");
  GetRecordSetW("time").initializeTags();

  Require<int>("cycle", Tags::DEFAULT, "cycle");
  Require<int>("cycle", Tags::NEXT, "cycle");
  GetRecordSetW("cycle").initializeTags();

  Require<double>("dt", Tags::DEFAULT, "coordinator");
  GetRecordW("dt", Tags::DEFAULT, "coordinator").set_initialized();

  Require<int>("position", Tags::DEFAULT, "position");
  GetRecordSetW("position").initializeTags();

  Teuchos::OSTab tab = vo_->getOSTab();

  // Ensure compatibility of all the evaluators -- each evaluator's dependencies
  // must provide what is required of that evaluator.
  for (auto& e : evaluators_) {
    for (auto& r : e.second) {
      if (!r.second->ProvidesKey(e.first, r.first)) {
        Errors::Message msg;
        msg << "Evaluator \"" << e.first << "\" with tag \"" << r.first.get()
            << "\" does not provide its own key.";
        Exceptions::amanzi_throw(msg);
      }

      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        Teuchos::OSTab tab1 = vo_->getOSTab();
        *vo_->os() << "checking compatibility: \"" << e.first << "\" + \"" << r.first.get() << "\"\n";
      }
      r.second->EnsureCompatibility(*this);
    }
  }

  // -- Create the data for all fields.
  for (auto& r : data_) {
    r.second->CreateData();

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "RecordSet \"" << r.first << "\", tags: ";

      auto& field = GetRecord(r.first);
      for(auto& e : *r.second) *vo_->os() << "\"" << e.first.get() << "\" ";
      if (field.ValidType<CompositeVector>()) {
        const auto& cv = Get<CompositeVector>(r.first);
        *vo_->os() << "comps: ";
        for (auto comp = cv.begin(); comp != cv.end(); ++comp) *vo_->os() << *comp << " ";
      } else {
        *vo_->os() << "not CV";
      }
      *vo_->os() << "\n";
    }
  }

  for (auto& deriv : derivs_) {
    // Some PKs allows an evalutator to be either independent or secondary.
    auto type = GetEvaluator(deriv.first).get_type();
    if (type == EvaluatorType::SECONDARY || type == EvaluatorType::PRIMARY) {
      deriv.second->CreateData();
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Derivative: D" << deriv.first << "/Dtags" << std::endl;
    }
  }

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab1 = vo_->getOSTab();
    *vo_->os() << "Setup is complete.\n\n";
  }
}


void State::Initialize()
{
  // Set metadata
  GetW<int>("cycle", Tags::DEFAULT, "cycle") = 0;
  data_.at("cycle")->initializeTags();

  GetW<double>("time", Tags::DEFAULT, "time") = 0.0;
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
      auto owner = GetRecord(e.first).owner();
      auto& field = GetRecordW(e.first, owner);
      const auto& copy = S->GetRecord(e.first);

      if (copy.initialized()) {
        if (field.ValidType<double>()) {
          field.Set<double>(owner, copy.Get<double>());
        } else if (field.ValidType<int>()) {
          field.Set<int>(owner, copy.Get<int>());
        } else if (field.ValidType<std::vector<double> >()) {
          field.Set<std::vector<double>>(owner, copy.Get<std::vector<double>>());
        } else if (field.ValidType<CompositeVector>()) {
          field.Set<CompositeVector>(owner, copy.Get<CompositeVector>());
        } else if (field.ValidType<AmanziGeometry::Point>()) {
          field.Set<AmanziGeometry::Point>(owner, copy.Get<AmanziGeometry::Point>());
        } else {
          std::stringstream ss;
          ss << "Trying to copy \"" << e.first << "\" with unknown type.\n";
          Errors::Message msg(ss.str());
          Exceptions::amanzi_throw(msg);
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
      *vo_->os() << "initializing eval: \"" << e.first << "\" with 1st tag \"" << tag->first.get() << "\"" << std::endl;
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

  Tag failed;
  if (state_plist_.isSublist("initial conditions")) {
    for (auto& e : data_) {
      bool flag(false);
      if (pre_initialization || !e.second->isInitialized(failed)) {
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
    ss << "Evaluator for field \"" << key << "\" at tag \"" << tag.get()
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
    ss << "Evaluator for field \"" << key << "\" at tag \"" << tag.get()
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

  return FEList().sublist(key);
}

}  // namespace Amanzi
