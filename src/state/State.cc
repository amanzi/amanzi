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

// Amanzi
#include "CompositeVector.hh"
#include "DomainSet.hh"
#include "errors.hh"
#include "Mesh.hh"
#include "StringExt.hh"

// Amanzi::State
#include "EvaluatorPrimary.hh"
#include "Evaluator_Factory.hh"
#include "StateDefs.hh"
#include "State.hh"
#include "Checkpoint.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructors
// -----------------------------------------------------------------------------
State::State()
  : state_plist_(Teuchos::rcp(new Teuchos::ParameterList("state"))),
    vo_(Teuchos::rcp(new VerboseObject("State", "none"))){};

State::State(const Teuchos::RCP<Teuchos::ParameterList>& state_plist)
  : state_plist_(state_plist), vo_(Teuchos::rcp(new VerboseObject("State", *state_plist)))
{
  // touch the "evaluators" list to make sure it exists, even in tests that
  // don't use it, to make sure that const calls to FEList() can work.
  state_plist_->sublist("evaluators");
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
    std::stringstream messagestream;
    messagestream << "Mesh " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return mesh;
};


Teuchos::RCP<AmanziMesh::Mesh>
State::GetDeformableMesh(Key key)
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


bool
State::IsDeformableMesh(const Key& key) const
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


Teuchos::RCP<AmanziMesh::Mesh>
State::GetMesh_(const Key& key) const
{
  if (key.empty()) return GetMesh_("domain");

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
State::RequireEvaluator(const Key& key, const Tag& tag)
{
  CheckIsDebugEval_(key, tag);

  // does it already exist?
  if (HasEvaluator(key, tag)) return GetEvaluator(key, tag);

  // See if the key is provided by another existing evaluator.
  for (auto& e : evaluators_) {
    if (Keys::hasKey(e.second, tag) && e.second.at(tag)->ProvidesKey(key, tag)) {
      auto& evaluator = e.second.at(tag);
      SetEvaluator(key, tag, evaluator);
      return *evaluator;
    }
  }

  // Check if this should be an aliased evaluator
  //
  // NOTE: rather than a pointer to the other evaluator, likely this should be
  // a special class that puts some guardrails on to make sure this evaluator
  // is ONLY used when the times match, or some data is valid, or some other
  // constraint.  This should not be used during this tag's subcycling, only
  // after subcycling/during syncronization/during another tag's
  // subcycling. --ETC
  if (tag == Tags::NEXT && evaluators_.count(key)) {
    for (const auto& other_tag : evaluators_.at(key)) {
      if (Keys::in(other_tag.first.get(), "next")) {
        // alias!
        SetEvaluator(key, tag, other_tag.second);
        return *other_tag.second;
      }
    }
  }

  // Create the evaluator from State's plist
  if (HasEvaluatorList(key)) {
    // -- Get this evaluator's plist.
    auto sublist = GetEvaluatorListPtr_(key);
    std::string evaluator_type = sublist->get<std::string>("evaluator type");
    bool by_region = Keys::ends_with(evaluator_type, "by region");
    bool initial_by_region = by_region;

    // -- move any model pars strings into additional model pars
    if (sublist->isParameter("model parameters") &&
        sublist->isType<std::string>("model parameters")) {
      if (sublist->isParameter("additional model parameters")) {
        if (sublist->isType<std::string>("additional model parameters")) {
          std::string add_par = sublist->get<std::string>("additional model parameters");
          sublist->remove("additional model parameters");
          sublist->set<Teuchos::Array<std::string>>("additional model parameters", std::vector<std::string>{add_par,});
        }

        auto add_pars = sublist->get<Teuchos::Array<std::string>>("additional model parameters");
        add_pars.push_back(sublist->get<std::string>("model parameters"));
        sublist->remove("additional model parameters");
        sublist->set<Teuchos::Array<std::string>>("additional model parameters", add_pars);
      } else {
        sublist->set<Teuchos::Array<std::string>>("additional model parameters",
                std::vector<std::string>{sublist->get<std::string>("model parameters")});
      }
      sublist->remove("model parameters");

    } else if (sublist->isParameter("model parameters") &&
               sublist->isType<Teuchos::Array<std::string>>("model parameters")) {
      if (sublist->isParameter("additional model parameters")) {
        if (sublist->isType<std::string>("additional model parameters")) {
          std::string add_par = sublist->get<std::string>("additional model parameters");
          sublist->remove("additional model parameters");
          sublist->set<Teuchos::Array<std::string>>("additional model parameters", std::vector<std::string>{add_par,});
        }

        auto add_pars = sublist->get<Teuchos::Array<std::string>>("additional model parameters");
        for (const auto& par : sublist->get<Teuchos::Array<std::string>>("model parameters")) {
          add_pars.push_back(par);
        }
        sublist->remove("additional model parameters");
        sublist->set<Teuchos::Array<std::string>>("additional model parameters", add_pars);
      } else {
        sublist->set<Teuchos::Array<std::string>>("additional model parameters",
                sublist->get<Teuchos::Array<std::string>>("model parameters"));
      }
      sublist->remove("model parameters");
    }

    // update model parameters with additional model parameters
    if (sublist->isParameter("additional model parameters")) {
      if (sublist->isType<std::string>("additional model parameters")) {
        std::string add_par = sublist->get<std::string>("additional model parameters");
        sublist->remove("additional model parameters");
        sublist->set<Teuchos::Array<std::string>>("additional model parameters", std::vector<std::string>{add_par,});
      }

      auto add_pars = sublist->get<Teuchos::Array<std::string>>("additional model parameters");
      for (const auto& add_par : add_pars) {
        sublist->sublist("model parameters").setParametersNotAlreadySet(GetModelParameters(add_par));
      }
      sublist->remove("additional model parameters");
    }

    if (sublist->isSublist("model parameters")) {
      Teuchos::ParameterList& model_pars = sublist->sublist("model parameters");

      // are we by region?
      for (const auto& p : model_pars) {
        if (model_pars.isSublist(p.first)) {
          by_region = true;
          break;
        }
      }

      if (by_region) {
        // explode the list, putting parameters in each region
        for (const auto& entry : model_pars) {
          if (!model_pars.isSublist(entry.first)) {
            for (auto& region_list : model_pars) {
              if (model_pars.isSublist(region_list.first)) {
                model_pars.sublist(region_list.first).setEntry(entry.first, entry.second);
              }
            }
            model_pars.remove(entry.first);
          }
        }

        Evaluator_Factory eval_fac;
        if (!initial_by_region &&
            eval_fac.HasEntry(evaluator_type+" by region")) {
          sublist->set<std::string>("evaluator type", evaluator_type+" by region");
        }
      }
    }

    // -- Create and set the evaluator.
    Evaluator_Factory evaluator_factory;
    sublist->set("tag", tag.get());
    auto evaluator = evaluator_factory.createEvaluator(sublist);
    SetEvaluator(key, tag, evaluator);
    return *evaluator;
  }

  // is it a meshed quantity?
  // recursive calls will result in the above, HasEvaluatorList() branch being taken.
  if (Keys::getVarName(key) == "cell_volume") {
    auto& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "cell volume");
    return RequireEvaluator(key, tag);
  } else if (Keys::getVarName(key) == "face_area") {
    auto& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "face area");
    return RequireEvaluator(key, tag);
  } else if (Keys::getVarName(key) == "elevation") {
    auto& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "elevation");
    return RequireEvaluator(key, tag);
  } else if (Keys::getVarName(key) == "slope_magnitude") {
    auto& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "slope magnitude");
    return RequireEvaluator(key, tag);
  } else if (Keys::getVarName(key) == "aspect") {
    auto& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "aspect");
    return RequireEvaluator(key, tag);
  } else if (Keys::getVarName(key) == "mesh") {
    auto& cv_list = GetEvaluatorList(key);
    cv_list.set("evaluator type", "static mesh");
    return RequireEvaluator(key, tag);
  }

  // cannot find the evaluator, error
  Errors::Message message;
  message << "Evaluator \"" << key << "@" << tag.get() << "\" cannot be created in State. "
          << "Verify (1) SetEvaluator is called or (2) name exists in state->evaluators.";
  Exceptions::amanzi_throw(message);
  return *Evaluator_Factory().createEvaluator(Teuchos::null); // silences warning
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


void
State::WriteDependencyGraph() const
{
  // // FIXME -- this is not what it used to be.  This simply writes data
  // // struture, and is not the dependency graph information at all.  Rename
  // // this, then recover the old WriteDependencyGraph method, which wrote a list
  // // of all evaluators and their dependnecies that could be read in networkx
  // // for plotting the dag. --ETC

  // if (vo_->os_OK(Teuchos::VERB_HIGH)) {
  //   *vo_->os() << "------------------------------------------" << std::endl
  //              << "Dependency & Structure list for evaluators" << std::endl
  //              << "------------------------------------------" << std::endl;
  //   for (auto& e : evaluators_) {
  //     for (auto& r : e.second) *vo_->os() << *r.second;
  //     if (GetRecord(e.first, e.second.begin()->first).ValidType<CompositeVector>()) {
  //       Teuchos::OSTab tab1 = vo_->getOSTab();
  //       Teuchos::OSTab tab2 = vo_->getOSTab();
  //       data_.at(e.first)->GetFactory<CompositeVector, CompositeVectorSpace>().print(*vo_->os());
  //       *vo_->os() << std::endl;
  //     }
  //   }

  //   *vo_->os() << "------------------------------" << std::endl
  //              << "Structure list for derivatives" << std::endl
  //              << "------------------------------" << std::endl;
  //   for (auto& e : derivs_) {
  //     *vo_->os() << "D" << e.first << "/D{ ";
  //     for (const auto& wrt : *e.second) *vo_->os() << wrt.first.get() << " ";
  //     *vo_->os() << "}" << std::endl;
  //     auto wrt_tag = e.second->begin();
  //     if (e.second->GetRecord(wrt_tag->first).ValidType<CompositeVector>()) {
  //       Teuchos::OSTab tab1 = vo_->getOSTab();
  //       Teuchos::OSTab tab2 = vo_->getOSTab();
  //       e.second->GetFactory<CompositeVector, CompositeVectorSpace>().print(*vo_->os());
  //       *vo_->os() << std::endl;
  //     }
  //   }
  // }
}


void
State::WriteStatistics(Teuchos::Ptr<const VerboseObject> vo,
                       const Teuchos::EVerbosityLevel vl) const
{
  if (vo == Teuchos::null) vo = vo_.ptr();
  Teuchos::OSTab tab = vo->getOSTab();

  // sort data in alphabetic order
  std::set<std::string> sorted;
  for (auto it = data_begin(); it != data_end(); ++it) { sorted.insert(it->first); }

  if (vo->os_OK(vl)) {
    *vo->os() << std::endl << "Field                                    Min/Max/Avg" << std::endl;

    for (auto name : sorted) { GetRecordSet(name).WriteStatistics(*vo); }
  }
}


// -----------------------------------------------------------------------------
// State handles model parameters.
// -----------------------------------------------------------------------------
const Teuchos::ParameterList&
State::GetModelParameters(const std::string& modelname)
{
  Teuchos::ParameterList& model_params = state_plist_->sublist("model parameters");
  if (!model_params.isSublist(modelname)) {
    Errors::Message message;
    message << "No set of \"model parameters\" named \"" << modelname
            << "\" in state->model parameters.";
    Exceptions::amanzi_throw(message);
  }
  return model_params.sublist(modelname);
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
  // NOTE: Need to add a check for loops as this DAG is created, otherwise the
  // user can create a loop in the input file that will seg-fault the
  // code. --ETC
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
        for (const auto& comp : *cvs) *vo_->os() << comp << " ";
      } else {
        *vo_->os() << "not CV";
      }
      *vo_->os() << "\n";
    }
  }

  // Create the data for all derivatives
  for (auto& deriv : derivs_) { deriv.second->CreateData(); }

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
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "initializing field copies..." << std::endl;
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
State::InitializeFields()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  // Initialize through constants list
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "initializing data through constants list..." << std::endl;
  }

  Teuchos::ParameterList& constants_list = state_plist_->sublist("constants");

  for (auto& e : data_) {
    std::string flag("... skipped");
    Tag not_init_tag;
    if (!e.second->isInitialized(not_init_tag)) {
      if (constants_list.isSublist(e.first)) {
        flag = "[ok]";
        e.second->Initialize(constants_list.sublist(e.first));
      } else {
        // check for domain set
        KeyTriple split;
        bool is_ds = Keys::splitDomainSet(e.first, split);
        Key ds_name = std::get<0>(split);
        if (is_ds) {
          Key lifted_key = Keys::getKey(ds_name, "*", std::get<2>(split));
          if (constants_list.isSublist(lifted_key)) {
            flag = "[ok]";
            // make a copy
            Teuchos::ParameterList sublist = constants_list.sublist(lifted_key);
            sublist.setName(e.first);
            e.second->Initialize(sublist);
          }
        }
      }

      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
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
      std::stringstream ss;
      ss << "Variable \"" << e.first << "\" with tag \"" << failed.get()
         << "\" was not initialized\n";
      Errors::Message msg(ss.str());
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
    state_plist_->get<Teuchos::Array<std::string>>("blacklist", empty).toVector();

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
    state_plist_->get<Teuchos::Array<std::string>>("whitelist", empty).toVector();

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
  } catch (std::out_of_range) {
    std::stringstream ss;
    ss << "Evaluator for field \"" << key << "\" at tag \"" << tag
       << "\" does not exist in the state.";
    Errors::Message message(ss.str());
    throw(message);
  }
}


Teuchos::RCP<Evaluator>
State::GetEvaluatorPtr(const Key& key, const Tag& tag)
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
void
State::SetEvaluator(const Key& key, const Tag& tag, const Teuchos::RCP<Evaluator>& evaluator)
{
  CheckIsDebugEval_(key, tag);
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
      if (FEList().isParameter(lifted_key)) { return FEList().sublist(lifted_key); }
    }
  }

  // return an empty new list
  return FEList().sublist(key);
}


bool
State::HasEvaluatorList(const Key& key) const
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


Teuchos::RCP<Teuchos::ParameterList>
State::GetEvaluatorListPtr_(const Key& key)
{
  if (FEList().isSublist(key)) {
    return Teuchos::sublist(Teuchos::sublist(state_plist_, "evaluators"), key);
  } else {
    KeyTriple split;
    bool is_ds = Keys::splitDomainSet(key, split);
    if (is_ds) {
      Key lifted_key = Keys::getKey(std::get<0>(split), "*", std::get<2>(split));
      if (FEList().isSublist(lifted_key)) {
        auto list_copy = Teuchos::rcp(new Teuchos::ParameterList(GetEvaluatorList(lifted_key)));
        list_copy->setName(key);
        // NOTE: should we put list_copy in "evaluators" list?
        return list_copy;
      }
    }
  }
  return Teuchos::null;
}


void
State::CheckIsDebugEval_(const Key& key, const Tag& tag)
{
  // check for debugging.  This provides a line for setting breakpoints for
  // debugging PK and Evaluator dependencies.
#ifdef ENABLE_DBC
  Teuchos::Array<Key> debug_evals = state_plist_->sublist("debug").get<Teuchos::Array<std::string>>(
    "evaluators", Teuchos::Array<Key>());
  if (std::find(debug_evals.begin(), debug_evals.end(), key) != debug_evals.end()) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "State: Evaluator for debug field \"" << key << "@" << tag << "\" was required."
                 << std::endl;
    }
    if (tag == Tags::DEFAULT) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) { *vo_->os() << " -- default tag" << std::endl; }
    }
  }
#endif
}

void
State::CheckIsDebugData_(const Key& key, const Tag& tag)
{
  // check for debugging.  This provides a line for setting breakpoints for
  // debugging PK and Evaluator dependencies.
#ifdef ENABLE_DBC
  Teuchos::Array<Key> debug_evals =
    state_plist_->sublist("debug").get<Teuchos::Array<std::string>>("data", Teuchos::Array<Key>());
  if (std::find(debug_evals.begin(), debug_evals.end(), key) != debug_evals.end()) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "State: data for debug field \"" << key << "@" << tag << "\" was required."
                 << std::endl;
    }
    if (tag == Tags::DEFAULT) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) { *vo_->os() << " -- default tag" << std::endl; }
    }
  }
#endif
}

} // namespace Amanzi
