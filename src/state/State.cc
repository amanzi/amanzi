/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

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
#include "EvaluatorAlias.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "Checkpoint.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------
State::State()
{
  vo_ = Teuchos::rcp(new VerboseObject("State", "none"));
};

State::State(Teuchos::ParameterList& state_plist)
  : state_plist_(state_plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("State", state_plist_));

  // touch the "evaluators" list to make sure it exists, even in tests that
  // don't use it, to make sure that const calls to FEList() can work.
  state_plist_.sublist("evaluators");
};


// -----------------------------------------------------------------------------
// State handles mesh management.
// -----------------------------------------------------------------------------
void
State::RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh, bool deformable)
{
  RegisterMesh("domain", mesh, deformable);
}


void
State::RegisterMesh(const Key& key, const Teuchos::RCP<AmanziMesh::Mesh>& mesh, bool deformable)
{
  meshes_.insert(std::make_pair(key, std::make_pair(mesh, deformable)));
};


void
State::AliasMesh(const Key& target, const Key& alias)
{
  bool deformable = IsDeformableMesh(target);
  Teuchos::RCP<AmanziMesh::Mesh> mesh = GetMesh_(target);
  RegisterMesh(alias, mesh, deformable);
  mesh_aliases_[alias] = target;

  if (GetMesh_(target + "_3d") != Teuchos::null) {
    RegisterMesh(alias + "_3d", GetMesh_(target + "_3d"), deformable);
    mesh_aliases_[alias + "_3d"] = target + "_3d";
  }
};


bool
State::IsAliasedMesh(const Key& key) const
{
  return Keys::hasKey(mesh_aliases_, key);
}


Teuchos::RCP<AmanziMesh::Mesh>
State::GetMesh(const Key& key) const
{
  Teuchos::RCP<AmanziMesh::Mesh> mesh;
  if (key.empty()) {
    mesh = GetMesh_("domain");
  } else {
    mesh = GetMesh_(key);
  }
  if (mesh == Teuchos::null) {
    Errors::Message msg;
    msg << "Mesh \"" << key << "\" does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  return mesh;
};


Teuchos::RCP<AmanziMesh::Mesh>
State::GetDeformableMesh(Key key)
{
  if (key.empty() ) key = "domain";

  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    if (lb->second.second) {
      return lb->second.first;
    } else {
      Errors::Message msg;
      msg << "Mesh " << key << " is not deformable.";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg;
    msg << "Mesh \"" << key << "\" does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  return Teuchos::null;
}


bool
State::IsDeformableMesh(const Key& key) const
{
  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second.second;
  } else {
    Errors::Message msg;
    msg << "Mesh \"" << key << "\" does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  return false;
}


Teuchos::RCP<AmanziMesh::Mesh>
State::GetMesh_(const Key& key) const
{
  if (key.empty() ) return GetMesh_("domain");

  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second.first;
  } else {
    return Teuchos::null;
  }
}


void
State::RegisterDomainSet(const Key& name, const Teuchos::RCP<AmanziMesh::DomainSet> set)
{
  domain_sets_[name] = set;
}

bool
State::HasDomainSet(const Key& name) const
{
  return (bool)domain_sets_.count(name);
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
// State handles data
// -----------------------------------------------------------------------------
//
// RecordSets
//
bool
State::HasRecordSet(const Key& fieldname) const
{
  return Keys::hasKey(data_, fieldname);
}

const RecordSet&
State::GetRecordSet(const Key& fieldname) const
{
  if (!HasRecordSet(fieldname)) {
    Errors::Message msg;
    msg << "State does not have a RecordSet named \"" << fieldname << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return *data_.at(fieldname);
}

RecordSet&
State::GetRecordSetW(const Key& fieldname)
{
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
bool
State::HasRecord(const Key& fieldname, const Tag& tag) const
{
  return HasRecordSet(fieldname) && GetRecordSet(fieldname).HasRecord(tag);
}

const Record&
State::GetRecord(const Key& fieldname, const Tag& tag) const
{
  return GetRecordSet(fieldname).GetRecord(tag);
}

Record&
State::GetRecordW(const Key& fieldname, const Key& owner)
{
  auto& r = GetRecordSetW(fieldname).GetRecord(Tags::DEFAULT);
  r.AssertOwnerOrDie(owner);
  return r;
}
Record&
State::GetRecordW(const Key& fieldname, const Tag& tag, const Key& owner)
{
  auto& r = GetRecordSetW(fieldname).GetRecord(tag);
  r.AssertOwnerOrDie(owner);
  return r;
}

//
// Derivative RecordSets
//
bool
State::HasDerivativeSet(const Key& fieldname, const Tag& tag) const
{
  return Keys::hasKey(derivs_, Keys::getKey(fieldname, tag));
}

const RecordSet&
State::GetDerivativeSet(const Key& fieldname, const Tag& tag) const
{
  if (!HasDerivativeSet(fieldname, tag)) {
    Errors::Message msg;
    msg << "State does not have a Derivative RecordSet for field \"" << fieldname << "\" at tag \""
        << tag.get() << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return *derivs_.at(Keys::getKey(fieldname, tag));
}

RecordSet&
State::GetDerivativeSetW(const Key& fieldname, const Tag& tag)
{
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
bool
State::HasDerivative(const Key& key, const Tag& tag, const Key& wrt_key, const Tag& wrt_tag) const
{
  auto keytag = Keys::getKey(key, tag);
  if (Keys::hasKey(derivs_, keytag)) {
    Tag der_tag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    return derivs_.at(keytag)->HasRecord(der_tag);
  }
  return false;
}
bool
State::HasDerivative(const Key& key, const Key& wrt_key) const
{
  return HasDerivative(key, Tags::DEFAULT, wrt_key, Tags::DEFAULT);
}


// -----------------------------------------------------------------------------
// State handles data evaluation.
// -----------------------------------------------------------------------------
Evaluator&
State::RequireEvaluator(const Key& key, const Tag& tag, bool alias_ok)
{
  bool debug = CheckIsDebugEval_(key, tag);
  Teuchos::OSTab tab = vo_->getOSTab();

  // does it already exist?
  if (HasEvaluator(key, tag)) {
    if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "  ... evaluator exists" << std::endl;
    }
    return GetEvaluator(key, tag);
  }

  // See if the key is provided by another existing evaluator.
  for (auto& e : evaluators_) {
    if (Keys::hasKey(e.second, tag) && e.second.at(tag)->ProvidesKey(key, tag)) {
      if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << "  ... evaluator is provided by evaluator for \""
                   << Keys::getKey(e.first, tag) << "\"" << std::endl;
      }
      auto& evaluator = e.second.at(tag);
      SetEvaluator_(key, tag, evaluator);
      evaluator->EnsureEvaluators(*this);
      return *evaluator;
    }
  }

  if (alias_ok) {
    // Check if this should be an aliased evaluator
    //
    // NOTE: This branch gets used, for example, when an evaluator is created
    // to do an observation, which uses (or is) a variable only computed at a
    // subcycled tag.
    //
    if (tag == Tags::NEXT && evaluators_.count(key)) {
      std::vector<Tag> other_nexts;
      for (const auto& other_tag : evaluators_.at(key)) {
        // can only auto-alias when there is only one other next tag
        if (Keys::in(other_tag.first.get(), "next")) {
          // BROKEN -- don't alias to something that depends upon this! --ETC
          other_nexts.push_back(other_tag.first);
        }
      }
      if (other_nexts.size() == 1) {
        Tag other_next = other_nexts[0];
        if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          *vo_->os() << "  ... evaluator aliased to \"" << other_next << "\"" << std::endl;
        }
        // alias!  Note the EvaluatorAlias makes sure that the data is shared
        // as well as the evaluator during its EnsureCompatibility()
        auto evaluator = Teuchos::rcp(new EvaluatorAlias(key, tag, other_next));
        SetEvaluator_(key, tag, evaluator);
        return *evaluator;
      }
    }
  }

  // Create the evaluator from State's plist
  // -- Get the Field Evaluator plist
  Teuchos::ParameterList& fm_plist = state_plist_.sublist("evaluators");
  Teuchos::ParameterList fe_plist;
  bool found = false;
  if (HasEvaluatorList(Keys::getKey(key, tag, true))) {
    // -- evaluator for just this key@tag combination
    fe_plist = GetEvaluatorList(Keys::getKey(key, tag));
    found = true;
    if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "  ... evaluator being created by tag-specific list" << std::endl;
    }
  } else if (HasEvaluatorList(key)) {
    // -- general evaluator for all tags
    fe_plist = GetEvaluatorList(key);
    found = true;
    if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "  ... evaluator being created by tag-generic list" << std::endl;
    }
  }
  if (found) {
    fe_plist.setName(key);
    // -- Check for special case -- alias
    if (fe_plist.isParameter("evaluator type") && fe_plist.isType<std::string>("evaluator type") &&
        fe_plist.get<std::string>("evaluator type") == "alias") {
      KeyTag target = Keys::splitKeyTag(fe_plist.get<std::string>("target"));
      if (target.first != key) {
        Errors::Message msg;
        msg << "Aliased evaluators may only alias across tags, not across keys.  Evaluator for \""
            << key << "\" does not match target key \"" << target.first << "\"";
        Exceptions::amanzi_throw(msg);
      }

      auto evaluator = Teuchos::rcp(new EvaluatorAlias(key, tag, target.second));
      SetEvaluator_(key, tag, evaluator);
      return *evaluator;
    }

    // -- Insert any model parameters.
    if (fe_plist.isParameter("model parameters") &&
        fe_plist.isType<std::string>("model parameters")) {
      std::string modelname = fe_plist.get<std::string>("model parameters");
      Teuchos::ParameterList modellist = GetModelParameters(modelname);
      fe_plist.set(modelname, modellist);
    }

    // -- Create and set the evaluator.
    Evaluator_Factory evaluator_factory;
    fe_plist.set("tag", tag.get());
    auto evaluator = evaluator_factory.createEvaluator(fe_plist);
    SetEvaluator_(key, tag, evaluator);
    evaluator->EnsureEvaluators(*this);
    return *evaluator;
  }

  // is it cell volume?
  if (Keys::getVarName(key) == "cell_volume") {
    Teuchos::ParameterList& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "cell volume");
    // recursive call will result in the above, HasEvaluatorList() branch being taken.
    if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "  ... evaluator created as cell volume" << std::endl;
    }
    return RequireEvaluator(key, tag);
  } else if (Keys::getVarName(key) == "slope_magnitude" || Keys::getVarName(key) == "elevation") {
    Teuchos::ParameterList& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "meshed elevation");
    // recursive call will result in the above, HasEvaluatorList() branch being taken.
    if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "  ... evaluator created as cell volume" << std::endl;
    }
    return RequireEvaluator(key, tag);
  }

  // Amanzi prefers to put constant-in-time values and functions in the
  // "initial conditions" list... can we make an IndependentVariable or
  // IndependentVariableConstant from stuff in that list?
  {
    if (HasICList(key)) {
      const Teuchos::ParameterList& ic_list = GetICList(key);
      if (ic_list.isSublist("function")) {
        Teuchos::ParameterList& e_list = GetEvaluatorList(key);
        e_list.set<std::string>("evaluator type", "independent variable");
        e_list.set<bool>("constant in time", true);
        e_list.sublist("function") = ic_list.sublist("function");

        if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          *vo_->os() << "  ... (from ICs list)";
        }

        return RequireEvaluator(key, tag);
      } else if (ic_list.isType<double>("value")) {
        Teuchos::ParameterList& e_list = GetEvaluatorList(key);
        e_list.set<std::string>("evaluator type", "independent variable constant");
        e_list.set<double>("value", ic_list.get<double>("value"));

        if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
          *vo_->os() << "  ... (from ICs list)";
        }

        return RequireEvaluator(key, tag);
      }
    }
  }

  // cannot find the evaluator, error
  Errors::Message msg;
  msg << "Evaluator \"" << key << "\" @ \"" << tag.get() << "\" cannot be created in State. "
      << "Verify (1) SetEvaluator is called or (2) name exists in state->evaluators.";
  Exceptions::amanzi_throw(msg);
  return *Evaluator_Factory().createEvaluator(fm_plist); // silences warning
}


bool
State::HasEvaluator(const Key& key, const Tag& tag) const
{
  if (Keys::hasKey(evaluators_, key)) {
    return Keys::hasKey(evaluators_.at(key), tag);
  } else {
    return false;
  }
}


Teuchos::RCP<const Functions::MeshPartition>
State::GetMeshPartition(Key key)
{
  Teuchos::RCP<const Functions::MeshPartition> mp = GetMeshPartition_(key);
  if (mp == Teuchos::null) {
    Errors::Message msg;
    msg << "Mesh partition \"" << key << "\" does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  return mp;
}


Teuchos::RCP<const Functions::MeshPartition>
State::GetMeshPartition_(Key key)
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
          part_list.get<Teuchos::Array<std::string>>("region list");
        Teuchos::RCP<Functions::MeshPartition> mp =
          Teuchos::rcp(new Functions::MeshPartition(AmanziMesh::CELL, region_list.toVector()));
        mp->Initialize(GetMesh(mesh_name), -1);
        mp->Verify();
        mesh_partitions_[key] = mp;
        return mp;
      }
    }
    return Teuchos::null;
  }
}


void
State::require_time(const Tag& tag, const Key& owner)
{
  Require<double>("time", tag, owner, false);
  if (!HasEvaluator("time", tag)) {
    Teuchos::ParameterList time_plist = GetEvaluatorList("time");
    time_plist.set<std::string>("tag", tag.get());
    auto evaluator = Teuchos::rcp(new EvaluatorPrimary<double>(time_plist));
    SetEvaluator("time", tag, evaluator);
  }
}

void
State::set_time(const Tag& tag, double value)
{
  Assign("time", tag, "time", value);
  auto time_eval = GetEvaluatorPtr("time", tag);
  auto time_eval_pv = Teuchos::rcp_dynamic_cast<EvaluatorPrimary<double>>(time_eval);
  AMANZI_ASSERT(time_eval_pv != Teuchos::null);
  time_eval_pv->SetChanged();
}


void
State::WriteDependencyGraph() const
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
      }
      *vo_->os() << std::endl;
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
const Teuchos::ParameterList&
State::GetModelParameters(const std::string& modelname) const
{
  AMANZI_ASSERT(state_plist_.isSublist("model parameters"));
  const Teuchos::ParameterList& model_plist = state_plist_.sublist("model parameters");
  AMANZI_ASSERT(model_plist.isSublist(modelname));
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
  require_time(Tags::DEFAULT);
  GetRecordSetW("time").initializeTags();
  require_time(Tags::CURRENT);
  GetRecordSetW("time").initializeTags();

  require_cycle(Tags::DEFAULT);
  GetRecordSetW("cycle").initializeTags();

  Require<double>("dt", Tags::DEFAULT, "dt", false);
  GetRecordW("dt", Tags::DEFAULT, "dt").set_initialized();

  Require<int>("position", Tags::DEFAULT, "position", false);
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
  //
  // Note this is done only because EnsureEvaluators() is called only in
  // RequireEvaluators(), and not in SetEvaluators().  See #654
  { // scope for copy
    EvaluatorMap evaluators_copy(evaluators_);
    for (auto& e : evaluators_copy) {
      for (auto& r : e.second) {
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          Teuchos::OSTab tab1 = vo_->getOSTab();
          *vo_->os() << "verify evaluator (DAG): \"" << e.first << "\" @ \"" << r.first << "\""
                     << std::endl;
        }
        r.second->EnsureEvaluators(*this);
      }
    }
  }

  // Second pass calls EnsureCompatibility, which checks data consistency.
  // This pass does not modify the graph.
  for (auto& e : evaluators_) {
    for (auto& r : e.second) {
      if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
        Teuchos::OSTab tab1 = vo_->getOSTab();
        *vo_->os() << "verify evaluator (compatibility): \"" << e.first << "\" @ \"" << r.first
                   << "\"" << std::endl;
      }
      r.second->EnsureCompatibility(*this);
    }
  }

  // Create the data for all fields.
  for (auto& r : data_) {
    r.second->CreateData();

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "RecordSet \"" << r.first << "\", tags: ";

      for (auto& e : *r.second) *vo_->os() << "\"" << e.first.get() << "\" ";
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


void
State::Initialize()
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
void
State::Initialize(const State& other)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab1 = vo_->getOSTab();
    *vo_->os() << "copying fields to new state..." << std::endl;
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


void
State::InitializeEvaluators()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  for (const auto& e : evaluators_) {
    for (const auto& tag : e.second) {
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "initializing eval: \"" << e.first << "\" @ \"" << tag.first << "\""
                   << std::endl;
      }

      tag.second->Update(*this, "state");
    }
    GetRecordSetW(e.first).initializeTags();
  }
}


void
State::InitializeFieldCopies(const Tag& reference_tag)
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


void
State::InitializeFields(const Tag& tag)
{
  bool pre_initialization = false;
  VerboseObject vo("State", state_plist_);

  // this directly overlaps with "initial conditions" --> "varname" -->
  // "restart file" capability, but is done globally instead of by variable.
  // Can the be refactored to use that instead?
  if (state_plist_.isParameter("initialization filename")) {
    pre_initialization = true;
    std::string filename = state_plist_.get<std::string>("initialization filename");
    Amanzi::Checkpoint file_input(filename, *this);

    for (auto it = data_.begin(); it != data_.end(); ++it) {
      auto owner = GetRecord(it->first, tag).owner();
      auto& r = GetRecordW(it->first, tag, owner);
      if (r.ValidType<CompositeVector>()) {
        r.ReadCheckpoint(file_input, tag, it->second->subfieldnames());

        // this is pretty hacky -- why are these ICs not in the PK's list?  And
        // if they aren't owned by a PK, they should be independent variables
        // not primary variables.
        if (HasEvaluator(it->first, tag)) {
          Evaluator& fm = GetEvaluator(it->first, tag);
          auto tmp = dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>*>(&fm);
          if (tmp != nullptr) {
            tmp->SetChanged();
          }
        }
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
    double t_ini = state_plist_.sublist("initial conditions").get<double>("time", 0.0);
    for (auto& e : data_) {
      std::string flag("... skipped");
      if (pre_initialization || !e.second->isInitialized(failed)) {
        if (state_plist_.sublist("initial conditions").isSublist(e.first)) {
          flag = "[ok]";
          if (state_plist_.sublist("initial conditions").isSublist(e.first)) {
            Teuchos::ParameterList sublist =
              state_plist_.sublist("initial conditions").sublist(e.first);
            sublist.set<double>("time", t_ini);
            e.second->Initialize(sublist, true);
          }
        } else {
          // check for domain set
          KeyTriple split;
          bool is_ds = Keys::splitDomainSet(e.first, split);
          Key ds_name = std::get<0>(split);
          if (is_ds) {
            Key lifted_key = Keys::getKey(ds_name, "*", std::get<2>(split));
            if (state_plist_.sublist("initial conditions").isSublist(lifted_key)) {
              flag = "[ok]";
              Teuchos::ParameterList sublist =
                state_plist_.sublist("initial conditions").sublist(lifted_key);
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
bool
State::CheckAllFieldsInitialized()
{
  Tag failed;
  for (auto& e : data_) {
    if (!e.second->isInitialized(failed)) {
      Errors::Message msg;
      msg << "Variable \"" << e.first << "\" with tag \"" << failed.get()
          << "\" was not initialized\n";
      Exceptions::amanzi_throw(msg);
      return false;
    }
  }
  return true;
}


// Utility for setting vis flags
void
State::InitializeIOFlags()
{
  Teuchos::Array<std::string> empty;

  // removing fields from vis dump
  std::vector<std::string> blacklist =
    state_plist_.get<Teuchos::Array<std::string>>("blacklist", empty).toVector();

  for (auto it = data_begin(); it != data_end(); ++it) {
    bool io_block(false);
    for (int m = 0; m < blacklist.size(); ++m) {
      std::regex pattern(blacklist[m]);
      io_block |= std::regex_match(it->first, pattern);
    }
    if (io_block) {
      for (auto& r : *it->second) r.second->set_io_vis(!io_block);
    }
  }

  // adding fields to vis dump
  std::vector<std::string> whitelist =
    state_plist_.get<Teuchos::Array<std::string>>("whitelist", empty).toVector();

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


Evaluator&
State::GetEvaluator(const Key& key, const Tag& tag)
{
  return *GetEvaluatorPtr(key, tag);
}


const Evaluator&
State::GetEvaluator(const Key& key, const Tag& tag) const
{
  try {
    return *evaluators_.at(key).at(tag);
  } catch (const std::out_of_range&) {
    Errors::Message msg;
    msg << "Evaluator for field \"" << key << "\" at tag \"" << tag
        << "\" does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  // silence warnings
  return *evaluators_.at(key).at(tag);
}


Teuchos::RCP<Evaluator>
State::GetEvaluatorPtr(const Key& key, const Tag& tag)
{
  try {
    return evaluators_.at(key).at(tag);
  } catch (const std::out_of_range&) {
    Errors::Message msg;
    msg << "Evaluator for field \"" << key << "\" at tag \"" << tag
        << "\" does not exist in the state.";
    Exceptions::amanzi_throw(msg);
  }
  // silence warnings
  return evaluators_.at(key).at(tag);
}


// Key defines field, Tag defines particular copy of field
void
State::SetEvaluator(const Key& key, const Tag& tag, const Teuchos::RCP<Evaluator>& evaluator)
{
  bool debug = CheckIsDebugEval_(key, tag);
  if (debug && vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "  ... evaluator manually set" << std::endl;
  }
  SetEvaluator_(key, tag, evaluator);
}


void
State::SetEvaluator_(const Key& key, const Tag& tag, const Teuchos::RCP<Evaluator>& evaluator)
{
  evaluators_[key][tag] = evaluator;
  // NOTE: specifically NOT doing this here, because Amanzi tends to construct
  // its DAG manually, evaluator by evaluator, in a hard-coded way.  Calling
  // this here would mean that Amanzi's PKs would have to always build their
  // DAG from the bottom up, which is not currently the case.  Therefore
  // calling EnsureEvaluators() while setting evaluators (which ensures that
  // the sub-graph rooted at this evaluator is always valid) is only done when
  // calling RequireEvaluator().
  //   evaluator->EnsureEvaluators(*this);
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
State::HasEvaluatorList(const Key& key) const
{
  if (FEList() .isSublist(key)) return true;
  // check for domain set
  KeyTriple split;
  bool is_ds = Keys::splitDomainSet(key, split);
  if (is_ds) {
    Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
    if (FEList() .isSublist(lifted_key)) return true;
  }
  return false;
}


Teuchos::ParameterList&
State::GetICList(const Key& key)
{
  if (ICList().isParameter(key)) {
    return ICList().sublist(key);
  } else {
    // check for domain set
    KeyTriple split;
    bool is_ds = Keys::splitDomainSet(key, split);
    if (is_ds) {
      Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
      if (ICList().isParameter(lifted_key)) {
        return ICList().sublist(lifted_key);
      }
    }
  }

  // return an empty new list
  return ICList().sublist(key);
}


bool
State::HasICList(const Key& key) const
{
  if (ICList() .isSublist(key)) return true;
  // check for domain set
  KeyTriple split;
  bool is_ds = Keys::splitDomainSet(key, split);
  if (is_ds) {
    Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
    if (ICList() .isSublist(lifted_key)) return true;
  }
  return false;
}


bool
State::CheckIsDebugEval_(const Key& key, const Tag& tag)
{
  // check for debugging.  This provides a line for setting breakpoints for
  // debugging PK and Evaluator dependencies.
#ifdef ENABLE_DBC
  if (!state_plist_.isSublist("debug") ) return false;

  Teuchos::Array<Key> debug_evals = state_plist_.sublist("debug").get<Teuchos::Array<std::string>>(
    "evaluators", Teuchos::Array<Key>());
  if (std::find(debug_evals.begin(), debug_evals.end(), key) != debug_evals.end()) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "State: Evaluator for debug field \"" << key << "@" << tag << "\" was required."
                 << std::endl;
    }
    if (tag == Tags::DEFAULT) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << " -- default tag" << std::endl;
      }
    }
    return true;
  }
#endif
  return false;
}

bool
State::CheckIsDebugData_(const Key& key, const Tag& tag)
{
  // check for debugging.  This provides a line for setting breakpoints for
  // debugging PK and Evaluator dependencies.
#ifdef ENABLE_DBC
  if (!state_plist_.isSublist("debug") ) return false;

  Teuchos::Array<Key> debug_evals =
    state_plist_.sublist("debug").get<Teuchos::Array<std::string>>("data", Teuchos::Array<Key>());
  if (std::find(debug_evals.begin(), debug_evals.end(), key) != debug_evals.end()) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "State: data for debug field \"" << key << "@" << tag << "\" was required."
                 << std::endl;
    }
    if (tag == Tags::DEFAULT) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        *vo_->os() << " -- default tag" << std::endl;
      }
    }
    return true;
  }
#endif
  return false;
}

} // namespace Amanzi
