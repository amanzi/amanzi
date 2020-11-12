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
#include <map>
#include <ostream>
#include <regex>

#include "boost/algorithm/string/predicate.hpp"

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "Mesh.hh"
#include "DomainSet.hh"
#include "MeshPartition.hh"
#include "CompositeVector.hh"
#include "FieldEvaluator_Factory.hh"
#include "cell_volume_evaluator.hh"
#include "rank_evaluator.hh"

#include "State.hh"

namespace Amanzi {

State::State() {};

State::State(Teuchos::ParameterList& state_plist) :
    state_plist_(state_plist),
    time_(0.0),
    cycle_(0),
    position_in_tp_(TIME_PERIOD_START)
{};


// copy constructor:
// Create a new State with different data but the same values.
//
// Could get a better implementation with a CopyMode, see TransportState in
// Amanzi as an example.  I'm not sure its needed at this point, however.
State::State(const State& other, StateConstructMode mode) :
    state_plist_(other.state_plist_),
    meshes_(other.meshes_),
    field_factories_(other.field_factories_),
    time_(other.time_),
    position_in_tp_(other.position_in_tp_),
    cycle_(other.cycle_) {

  if (mode == STATE_CONSTRUCT_MODE_COPY_DATA) {
    for (FieldMap::const_iterator f_it=other.fields_.begin();
         f_it!=other.fields_.end(); ++f_it) {
      fields_[f_it->first] = f_it->second->Clone();
    }

    for (const auto& fm_it : other.field_evaluators_) {
      bool copied = false;
      for (auto& fm_my : field_evaluators_) {
        if (fm_my.second->ProvidesKey(fm_it.first)) {
          field_evaluators_[fm_it.first] = fm_my.second;
          copied = true;
          break;
        }
      }
      if (!copied)
        field_evaluators_[fm_it.first] = fm_it.second->Clone();
    }

  } else {
    for (FieldMap::const_iterator f_it=other.fields_.begin();
         f_it!=other.fields_.end(); ++f_it) {
      fields_[f_it->first] = f_it->second;
    }

    for (FieldEvaluatorMap::const_iterator fm_it=other.field_evaluators_.begin();
         fm_it!=other.field_evaluators_.end(); ++fm_it) {
      field_evaluators_[fm_it->first] = fm_it->second;
    }
  }

  // pointer copy MeshPartitions -- const after initing.
  for (MeshPartitionMap::const_iterator mp_it=other.mesh_partitions_.begin();
       mp_it!=other.mesh_partitions_.end(); ++mp_it) {
    mesh_partitions_[mp_it->first] = mp_it->second;
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


//  Assign a state's data from another state.  Note this
// implementation requires the State being copied has the same structure (in
// terms of fields, order of fields, etc) as *this.  This really means that
// it should be a previously-copy-constructed version of the State.  One and
// only one State should be instantiated and populated -- all other States
// should be copy-constructed from that initial State.
void State::AssignDomain(const State& other, const std::string& domain)
{
  if (this != &other) {
    for (auto f_it : other.fields_) {
      if (Keys::getDomain(f_it.first) == domain) {
        auto myfield = GetField_(f_it.first);
        auto otherfield = f_it.second;

        if (myfield->type() == COMPOSITE_VECTOR_FIELD) {
          myfield->SetData(*otherfield->GetFieldData());
        } else if (myfield->type() == CONSTANT_VECTOR) {
          myfield->SetData(*otherfield->GetConstantVectorData());
        } else if (myfield->type() == CONSTANT_SCALAR) {
          myfield->SetData(*otherfield->GetScalarData());
        }
      }
    }

    for (auto fm_it : other.field_evaluators_) {
      if (Keys::getDomain(fm_it.first) == domain) {
        auto myfm = GetFieldEvaluator_(fm_it.first);
        *myfm = *fm_it.second;
      }
    }
  }
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
  mesh_aliases_[target] = alias;

  if (GetMesh_(target+"_3d") != Teuchos::null) {
    RegisterMesh(alias+"_3d", GetMesh_(target+"_3d"), deformable);
    mesh_aliases_[target+"_3d"] = alias+"_3d";
  }
};

bool State::IsAliasedMesh(const Key& key) const {
  return (bool) mesh_aliases_.count(key);
}


void State::RemoveMesh(const Key& key) {
  meshes_.erase(key);
};


Teuchos::RCP<const AmanziMesh::Mesh> State::GetMesh(const Key& key) const {
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


Teuchos::RCP<AmanziMesh::Mesh> State::GetDeformableMesh(Key key) {
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
};


bool State::IsDeformableMesh(const Key& key) const {
  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second.second;
  } else {
    std::stringstream messagestream;
    messagestream << "Mesh " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return false;
};


Teuchos::RCP<AmanziMesh::Mesh> State::GetMesh_(const Key& key) const {
  if (key.empty()) return GetMesh_("domain");
  
  mesh_iterator lb = meshes_.lower_bound(key);
  if (lb != meshes_.end() && !(meshes_.key_comp()(key, lb->first))) {
    return lb->second.first;
  } else {
    return Teuchos::null;
  }
};


void State::RegisterDomainSet(const Key& name,
        const Teuchos::RCP<AmanziMesh::DomainSet> set) {
  domain_sets_[name] = set;
}

bool State::HasDomainSet(const Key& name) const {
  return (bool) domain_sets_.count(name);
}

Teuchos::RCP<const AmanziMesh::DomainSet>
State::GetDomainSet(const Key& name) const {
  return domain_sets_.at(name);
}


// -----------------------------------------------------------------------------
// State handles data evaluation.
// -----------------------------------------------------------------------------
Teuchos::RCP<FieldEvaluator>
State::RequireFieldEvaluator(Key key) {
  Teuchos::RCP<FieldEvaluator> evaluator = GetFieldEvaluator_(key);

  // Get the evaluator from state's Plist
  if (evaluator == Teuchos::null) {
    // -- Get the Field Evaluator plist
    Teuchos::ParameterList& fm_plist = FEList();
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
      FieldEvaluator_Factory evaluator_factory;
      evaluator = evaluator_factory.createFieldEvaluator(sublist);
      SetFieldEvaluator(key, evaluator);
    }
  }

  // check to see if we have a flyweight evaluator
  if (evaluator == Teuchos::null) {
    KeyTriple split;
    bool is_ds = Keys::splitDomainSet(key, split);
    if (is_ds) {
      Key ds_name = std::get<0>(split);
      if (HasDomainSet(std::get<0>(split))) {
        auto ds = GetDomainSet(ds_name);
        // The name is a domain set prefixed name, and we have a domain set
        // of that name.  Grab the parameter list for the set's list and use
        // that to construct the evaluator.
        // -- Get this evaluator's plist.
        Key lifted_key = Keys::getKey(ds_name+"_*",std::get<2>(split));

        Teuchos::ParameterList& fm_plist = FEList();
        if (fm_plist.isSublist(lifted_key)) {
          Teuchos::ParameterList sublist = fm_plist.sublist(lifted_key);
          sublist.set("evaluator name", key);
          sublist.setName(key);

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
          fm_plist.set(key, sublist);
          FieldEvaluator_Factory evaluator_factory;
          evaluator = evaluator_factory.createFieldEvaluator(sublist);
            
          SetFieldEvaluator(key, evaluator);
        }            
      } else {
        Errors::Message msg;
        msg << "Model for field \"" << key << "\" is on a DomainSet, but no DomainSet of name \"" << ds_name << "\" exists in State.";
        Exceptions::amanzi_throw(msg);
      }
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

  // cannot find the evaluator, error
  if (evaluator == Teuchos::null) {
    std::stringstream msg;
    msg << "\nModel for field \"" << key << "\" cannot be created in State.\n";
    // for (auto fe = field_evaluator_begin(); fe != field_evaluator_end(); ++fe) {
    //   msg << fe->first << ":\n" << fe->second;
    // }
    Errors::Message message(msg.str());
    Exceptions::amanzi_throw(message);
  }
  return evaluator;
}


Teuchos::RCP<FieldEvaluator>
State::RequireFieldEvaluator(Key key, Teuchos::ParameterList& plist) {
  Teuchos::RCP<FieldEvaluator> evaluator = GetFieldEvaluator_(key);

  // Create a new evaluator.
  if (evaluator == Teuchos::null) {
    // -- Create and set the evaluator.
    FieldEvaluator_Factory evaluator_factory;
    evaluator = evaluator_factory.createFieldEvaluator(plist);
    SetFieldEvaluator(key, evaluator);
  }
  return evaluator;
}


Teuchos::RCP<FieldEvaluator> State::GetFieldEvaluator(const Key& key) {
  Teuchos::RCP<FieldEvaluator> evaluator = GetFieldEvaluator_(key);
  if (evaluator == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Model for field " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return evaluator;
};


Teuchos::RCP<FieldEvaluator> State::GetFieldEvaluator_(const Key& key) {
  auto f_it = field_evaluators_.find(key);
  if (f_it != field_evaluators_.end()) {
    return f_it->second;;
  } else {
    // See if the key is provided by another existing evaluator.
    for (evaluator_iterator f_it = field_evaluator_begin();
         f_it != field_evaluator_end(); ++f_it) {
      if (f_it->second->ProvidesKey(key)) {
        SetFieldEvaluator(key, f_it->second);
        return f_it->second;
      }
    }
  }
  return Teuchos::null;
};


void State::SetFieldEvaluator(Key key, const Teuchos::RCP<FieldEvaluator>& evaluator) {
  AMANZI_ASSERT(field_evaluators_[key] == Teuchos::null);
  field_evaluators_[key] = evaluator;
};


Teuchos::RCP<const Functions::MeshPartition> State::GetMeshPartition(Key key) {
  Teuchos::RCP<const Functions::MeshPartition> mp = GetMeshPartition_(key);
  if (mp == Teuchos::null) {
    std::stringstream messagestream;
    messagestream << "Mesh partition " << key << " does not exist in the state.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }
  return mp;
};


Teuchos::RCP<const Functions::MeshPartition> State::GetMeshPartition_(Key key) {
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
};


void State::WriteDependencyGraph() const {
  auto vo = Teuchos::rcp(new VerboseObject("State", state_plist_)); 
  if (vo->os_OK(Teuchos::VERB_HIGH)) {
    std::ofstream os("dependency_graph.txt", std::ios::out);
    for (auto fe=field_evaluator_begin(); fe!=field_evaluator_end(); ++fe) {
      os << *fe->second;
    }
    os.close();
  }
}


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
    Teuchos::RCP<Field_Scalar> newfield = Teuchos::rcp(new Field_Scalar(fieldname, owner));
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
    field = Teuchos::rcp(new Field_ConstantVector(fieldname, Key("state"), dimension));
    fields_[fieldname] = field;
  } else {
    Teuchos::RCP<Field_ConstantVector> cv =
        Teuchos::rcp_static_cast<Field_ConstantVector>(field);
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
Teuchos::RCP<CompositeVectorSpace>
State::RequireField(Key fieldname, Key owner) {
  Teuchos::RCP<Field> field = CheckConsistent_or_die_(fieldname,
          COMPOSITE_VECTOR_FIELD, owner);

  if (field == Teuchos::null) {
    // Create the field and CV factory.
    field = Teuchos::rcp(new Field_CompositeVector(fieldname, owner));
    fields_[fieldname] = field;
    field_factories_[fieldname] = Teuchos::rcp(new CompositeVectorSpace());
  } else if (owner != Key("state")) {
    field->set_owner(owner);
  }

  return field_factories_[fieldname];
};

// Vector Field requires
Teuchos::RCP<CompositeVectorSpace>
State::RequireField(Key fieldname, Key owner,
                    const std::vector<std::vector<std::string> >& subfield_names) {
  Teuchos::RCP<Field> field = CheckConsistent_or_die_(fieldname,
          COMPOSITE_VECTOR_FIELD, owner);

  if (field == Teuchos::null) {
    // Create the field and CV factory.
    field = Teuchos::rcp(new Field_CompositeVector(fieldname, owner, subfield_names));
    fields_[fieldname] = field;
    field_factories_[fieldname] = Teuchos::rcp(new CompositeVectorSpace());
  } else if (owner != Key("state")) {
    field->set_owner(owner);
  }

  return field_factories_[fieldname];
};


void State::RequireTimeTag(Key timetag) {

  for (int i=0; i<copy_tag_.size(); ++i) {
    if (copy_tag_[i] == timetag) return;
  }
  copy_tag_.push_back(timetag);
}

bool State::HasTimeTag(Key timetag) { 

  for (int i=0; i<copy_tag_.size(); ++i) {
    if (copy_tag_[i] == timetag) return true;
  }
  return false;
}

bool State::HasFieldCopy(Key key, Key timetag) {

  Teuchos::RCP<const Field> field = GetField_(key);
  if (field == Teuchos::null) return false;
  if (HasTimeTag(timetag)) return field->HasCopy(timetag);
  else return false;
}

Teuchos::RCP<Field> State::GetFieldCopy(Key fieldname, Key timetag, Key pk_name) {

  Key field_owner = GetField(fieldname)->owner();
  Teuchos::RCP<Field> field = GetField(fieldname, field_owner);
  return field->GetCopy(timetag, pk_name);

}

Teuchos::RCP<const Field> State::GetFieldCopy(Key fieldname, Key timetag) const{

  Teuchos::RCP<const Field> field = GetField(fieldname);
  return field->GetCopy(timetag);

}

void State::SetFieldCopy(Key fieldname, Key timetag, Key pk_name, const Teuchos::RCP<Field>& fieldcopy) {

   Teuchos::RCP<Field> field = GetField(fieldname, pk_name);
   return field->SetCopy(timetag, fieldcopy);

}


void State::RequireFieldCopy(Key fieldname, Key timetag, Key copy_owner) {

  Key field_owner = GetField(fieldname)->owner();
  Teuchos::RCP<Field> field = GetField(fieldname, field_owner);
  field->RequireCopy(timetag, copy_owner);
}


void State::CopyField(Key fieldname, Key timetag, Key pk_name) {

  Key field_owner = GetField(fieldname)->owner();
  Teuchos::RCP<Field> field = GetField(fieldname, field_owner);
  Teuchos::RCP<Field> field_copy = field->GetCopy(timetag, pk_name);

  *field_copy->GetFieldData() = *field->GetFieldData();
}


Teuchos::RCP<const CompositeVector> State::GetFieldCopyData(Key fieldname, Key tag) const{
  if (tag=="default") 
    return GetField(fieldname)->GetFieldData();
  else
    return GetField(fieldname)->GetCopy(tag)->GetFieldData();
}


Teuchos::RCP<CompositeVector> State::GetFieldCopyData(Key fieldname, Key tag, Key pk_name) {

  if (tag=="default") {
    Teuchos::RCP<Field> field = GetField(fieldname, pk_name);
    return field->GetFieldData();
  }
  else {
    Key field_owner = GetField(fieldname)->owner();
    Teuchos::RCP<Field> field = GetField(fieldname, field_owner);
    return field->GetCopy(tag, pk_name)->GetFieldData();
  }
}


void State::RequireGravity() {
  int dim = 3;
  Key fieldname("gravity");
  RequireConstantVector(fieldname, dim);
  std::vector<Key> subfield_names(dim);
  subfield_names[0] = "x";
  if (dim > 1) subfield_names[1] = "y";
  if (dim > 2) subfield_names[2] = "z";
  Teuchos::RCP<Field> field = GetField_(fieldname);
  Teuchos::RCP<Field_ConstantVector> cvfield =
    Teuchos::rcp_dynamic_cast<Field_ConstantVector>(field, true);
  cvfield->set_subfield_names(subfield_names);
  cvfield->CreateData();
  cvfield->Initialize(state_plist_.sublist("initial conditions").sublist("gravity"));
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
    messagestream << "Field " << fieldname << " does not exist in the state, pk=" << pk_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  } else if (record->owner() != pk_name) {
    std::stringstream messagestream;
    messagestream << "PK \"" << pk_name
                  << "\" is attempting write access to field \"" << fieldname
                  << "\" which is owned by \"" << GetField_(fieldname)->owner() <<"\"";
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

void State::SetField(Key fieldname, Key pk_name,
                     const Teuchos::RCP<Field>& new_record) {
  Teuchos::RCP<const Field> old_record = GetField_(fieldname);
  if (old_record != Teuchos::null) {
    if (old_record->owner() != pk_name) {
      std::stringstream messagestream;
      messagestream << "PK " << pk_name
                    << " is attempting to overwrite field " << fieldname
                    << " which is owned by " << old_record->owner();
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
  }
  fields_[fieldname] = new_record;
}

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
                    const Teuchos::RCP<CompositeVector>& data) {
  GetField(fieldname, pk_name)->SetData(data);
};


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
void State::Setup() {
  auto vo = Teuchos::rcp(new VerboseObject("State", state_plist_)); 
  Teuchos::OSTab tab = vo->getOSTab();
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
      SetFieldEvaluator(rank_key, Teuchos::rcp(new RankEvaluator(plist)));
    }
  }

  // Ensure compatibility of all the evaluators -- each evaluator's dependencies must
  // provide what is required of that evaluator.
  for (FieldEvaluatorMap::iterator evaluator=field_evaluators_.begin();
       evaluator!=field_evaluators_.end(); ++evaluator) {
    if (!evaluator->second->ProvidesKey(evaluator->first)) {
      std::stringstream messagestream;
      messagestream << "Field Evaluator " << evaluator->first
                    << " does not provide its own key.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
    if (vo->os_OK(Teuchos::VERB_HIGH)) {
      Teuchos::OSTab tab1 = vo->getOSTab();
      *vo->os() << "checking compatibility: " << *evaluator->second;
    }
    evaluator->second->EnsureCompatibility(Teuchos::ptr(this));
  }

  // Create all data for vector fields.
  // -- First use factories to instantiate composite vectors.
  for (FieldFactoryMap::iterator fac_it=field_factories_.begin();
       fac_it!=field_factories_.end(); ++fac_it) {
    GetField_(fac_it->first)->SetData(Teuchos::rcp(new CompositeVector(*fac_it->second)));
  }

  // -- Now create the data for all fields.
  for (field_iterator f_it = field_begin(); f_it != field_end(); ++f_it) {
    if (!f_it->second->initialized()) {
      f_it->second->CreateData();
    }
    // std::cout<<"Field "<<f_it->first<<": ";
    // if (f_it->second->type() == Amanzi::COMPOSITE_VECTOR_FIELD) {
    //   auto com_vec = f_it->second->GetFieldData();
    //     for (CompositeVector::name_iterator comp=com_vec->begin();
    //          comp!=com_vec->end(); ++comp) std::cout<<*comp<<" ";
    //   std::cout<<"\n";
    // }
  }
};


void State::Initialize() {
  // Initialize any other fields from state plist.
  InitializeFields();

  // Ensure that non-evaluator-based fields are initialized.
  CheckNotEvaluatedFieldsInitialized();

  // Initialize other field evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  CheckAllFieldsInitialized();

  // Write dependency graph.
  WriteDependencyGraph();

  // Reset io_vis flags using blacklist and whitelist
  InitializeIOFlags_();
};


void State::Initialize(Teuchos::RCP<State> S) {
  auto vo = Teuchos::rcp(new VerboseObject("State", state_plist_)); 
  Teuchos::OSTab tab = vo->getOSTab();
  if (vo->os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab1 = vo->getOSTab();
    *vo->os() << "copying fields to new state.." << std::endl;
  }
  
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    Teuchos::RCP<Field> copy = S->GetField_(field->fieldname());

    if (copy != Teuchos::null) {
      if (field->type() != copy->type()) {
        std::stringstream messagestream;
        messagestream << "States has fields with the same name but different types\n";
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
      }
      if (copy->initialized()) {
        switch (field->type()) {
        case CONSTANT_SCALAR:
          *field->GetScalarData() = *copy->GetScalarData();
          break;
        case CONSTANT_VECTOR:
          *field->GetConstantVectorData() = *copy->GetConstantVectorData();
          break;
        case COMPOSITE_VECTOR_FIELD:
          *field->GetFieldData() = *copy->GetFieldData();
          break;
        default:
          Errors::Message message("Copy field with unknown type\n");
          Exceptions::amanzi_throw(message);
        }
        field->set_initialized();      
      }
    }   
  }

  // Initialize any other fields from state plist.
  InitializeFields();

  // Ensure that non-evaluator-based fields are initialized.
  // CheckNotEvaluatedFieldsInitialized();

  // Initialize other field evaluators.
  InitializeEvaluators();

  // Ensure everything is owned and initialized.
  // CheckAllFieldsInitialized();

  // Write dependency graph.
  WriteDependencyGraph();

  // Reset io_vis flags using blacklist and whitelist
  InitializeIOFlags_();
};


void State::InitializeEvaluators() {
  VerboseObject vo("State", state_plist_); 
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "initializing evaluators..." << std::endl;
  }
  for (evaluator_iterator f_it = field_evaluator_begin(); f_it != field_evaluator_end(); ++f_it) {
    f_it->second->HasFieldChanged(Teuchos::Ptr<State>(this), "state");
    fields_[f_it->first]->set_initialized();
  }
};


void State::InitializeFields() {
  bool pre_initialization = false;
  VerboseObject vo("State", state_plist_); 

  if (state_plist_.isParameter("initialization filename")) {
    pre_initialization = true;
    std::string filename = state_plist_.get<std::string>("initialization filename");
    Amanzi::HDF5_MPI file_input(GetMesh()->get_comm(), filename);
    file_input.open_h5file();
    for (FieldMap::iterator f_it = fields_.begin(); f_it != fields_.end(); ++f_it) {
      if (f_it->second->type() == COMPOSITE_VECTOR_FIELD) {
        bool read_complete = f_it->second->ReadCheckpoint(file_input);
        if (read_complete)
          f_it->second->set_initialized();
      }
    }
    file_input.close_h5file();
  }
   
  // Initialize through initial condition
  if (vo.os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "initializing fields through initial conditions..." << std::endl;
  }

  for (FieldMap::iterator f_it = fields_.begin(); f_it != fields_.end(); ++f_it) {
    if (pre_initialization || (!f_it->second->initialized())) {
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
bool State::CheckNotEvaluatedFieldsInitialized() {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    if (!HasFieldEvaluator(f_it->first)) {
      // first check and see if there is a FieldEvaluator, but we haven't yet used it.
      Teuchos::RCP<FieldEvaluator> found_eval;
      for (evaluator_iterator f_eval_it = field_evaluator_begin();
         f_eval_it != field_evaluator_end(); ++f_eval_it) {
        if (f_eval_it->second->ProvidesKey(f_it->first)) {
          found_eval = f_eval_it->second;
          break;
        }
      }

      if (found_eval != Teuchos::null) {
        // found an evaluator that provides this key, set it.
        SetFieldEvaluator(f_it->first, found_eval);
      } else if (!field->initialized()) {
        // No evaluator, not intialized... FAIL.
        std::stringstream messagestream;
        messagestream << "Field " << field->fieldname() << " was not initialized. Owner:" << field->owner();
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
        return false;
      }
    }
  }
  return true;
};

// Make sure all fields that are not evaluated by a FieldEvaluator are
// initialized.  Such fields may be used by an evaluator field but are not in
// the dependency tree due to poor design.
bool State::CheckNotEvaluatedFieldsInitialized(Teuchos::RCP<State> S) {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    if (!HasFieldEvaluator(f_it->first)) {
      // first check and see if there is a FieldEvaluator, but we haven't yet used it.
      Teuchos::RCP<FieldEvaluator> found_eval;
      for (evaluator_iterator f_eval_it = field_evaluator_begin();
         f_eval_it != field_evaluator_end(); ++f_eval_it) {
        if (f_eval_it->second->ProvidesKey(f_it->first)) {
          found_eval = f_eval_it->second;
          break;
        }
      }

      if (found_eval != Teuchos::null) {
        // found an evaluator that provides this key, set it.
        SetFieldEvaluator(f_it->first, found_eval);
      } else if (!field->initialized()) {
        // No evaluator, not intialized... FAIL.
        std::stringstream messagestream;
        messagestream << "Field " << field->fieldname() << " was not initialized. Owner:" << field->owner();
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
        return false;
      }
    }
  }
  return true;
};


// Make sure all fields have gotten their IC, either from State or the owning PK.
bool State::CheckAllFieldsInitialized() {
  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    if (!field->initialized()) {
      // field was not initialized
      std::stringstream messagestream;
      messagestream << "Field " << field->fieldname() << " was not initialized. Owner:" << field->owner();
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
      return false;
    }
  }
  return true;
};



// Make sure all fields have gotten their IC, either from State or the owning PK.
// if state S has the same filed the data is copied

bool State::CheckAllFieldsInitialized(Teuchos::RCP<State> S) {


  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    Teuchos::RCP<Field> copy = S->GetField_(field->fieldname());
    if (copy != Teuchos::null) {
      if (field->type() != copy->type()) {
        std::stringstream messagestream;
        messagestream << "States has fields with the same name but different types\n";
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
        return false;
      }
      switch (field->type()) {
      case CONSTANT_SCALAR:
        *field->GetScalarData() = *copy->GetScalarData();
        break;
      case CONSTANT_VECTOR:
        *field->GetConstantVectorData() = *copy->GetConstantVectorData();
        break;
      case COMPOSITE_VECTOR_FIELD:
        *field->GetFieldData() = *copy->GetFieldData();
        break;
      default:
        Errors::Message message("Copy field with unknown type\n");
        Exceptions::amanzi_throw(message);
        return false;
      }
      field->set_initialized();
    }   
  }


  for (FieldMap::iterator f_it = fields_.begin();
       f_it != fields_.end(); ++f_it) {
    Teuchos::RCP<Field> field = f_it->second;
    if (!field->initialized()) {
      // field was not initialized
      std::stringstream messagestream;
      messagestream << "Field " << field->fieldname() << " was not initialized. Owner:" << field->owner();
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
  AMANZI_ASSERT(fieldname != Key(""));

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


// Utility for setting vis flags
void State::InitializeIOFlags_() {
  Teuchos::Array<std::string> empty;

  // removing fields from vis dump
  std::vector<std::string> blacklist = 
      state_plist_.get<Teuchos::Array<std::string> >("blacklist", empty).toVector();

  for (State::field_iterator field=field_begin(); field!=field_end(); ++field) {
    bool io_block(false);
    for (int m = 0; m < blacklist.size(); ++m) {
      std::regex pattern(blacklist[m]);
      io_block |= std::regex_match(field->first, pattern);
    }
    field->second->set_io_vis(!io_block);
  }

  // adding fields to vis dump
  std::vector<std::string> whitelist =
      state_plist_.get<Teuchos::Array<std::string> >("whitelist", empty).toVector();

  for (State::field_iterator field=field_begin(); field!=field_end(); ++field) {
    bool io_allow(false);
    for (int m = 0; m < whitelist.size(); ++m) {
      std::regex pattern(whitelist[m]);
      io_allow |= std::regex_match(field->first, pattern);
    }
    if (io_allow) field->second->set_io_vis(true);
  }
}


// Non-member function for vis.
void WriteVis(Visualization& vis,
              State& S) {
  if (!vis.is_disabled()) {
    // Create the new time step
    vis.CreateTimestep(S.time(), S.cycle(), vis.get_tag());

    // Write all fields to the visualization file, the fields know if they
    // need to be written.
    for (State::field_iterator field=S.field_begin(); field!=S.field_end(); ++field) {
      field->second->WriteVis(vis);
    }
    
    vis.WriteRegions();
    vis.WritePartition();

    // Finalize i/o.
    vis.FinalizeTimestep();
  }
};


// Non-member function for checkpointing.
void WriteCheckpoint(Checkpoint& chk,
                     State& S,
                     double dt,
                     bool final,
                     Amanzi::ObservationData* obs_data) {
  if ( !chk.is_disabled() ) {
    chk.CreateFile(S.cycle());

    // create hard link to the final file
    if (final && S.GetMesh()->get_comm()->MyPID() == 0)
      chk.CreateFinalFile(S.cycle());

    for (State::field_iterator field=S.field_begin(); field!=S.field_end(); ++field) {
      field->second->WriteCheckpoint(chk);
    }

    chk.WriteAttributes(S.time(), dt, S.cycle(), S.position());
    chk.WriteObservations(obs_data);
    
    chk.Finalize();
  }
};


// Non-member function for checkpointing.
double ReadCheckpoint(const Comm_ptr_type& comm,
                      State& S,
                      std::string filename) {
  HDF5_MPI checkpoint(comm, filename);
  checkpoint.open_h5file();

  // load the attributes
  double time(0.);
  checkpoint.readAttrReal(time, "time");
  S.set_time(time);

  double dt(0.);
  checkpoint.readAttrReal(dt, "dt");

  int cycle(0);
  checkpoint.readAttrInt(cycle, "cycle");
  S.set_cycle(cycle);

  int pos(0);
  checkpoint.readAttrInt(pos, "position");
  S.set_position(pos);

  // load the number of processes and ensure they are the same -- otherwise
  // the below just gives crap.
  int rank(-1);
  checkpoint.readAttrInt(rank, "mpi_comm_world_rank");
  if (comm->NumProc() != rank) {
    std::stringstream messagestream;
    messagestream << "Requested checkpoint file " << filename << " was created on "
                  << rank << " processes, making it incompatible with this run on "
                  << comm->NumProc() << " process.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // load the data
  for (State::field_iterator field=S.field_begin(); field!=S.field_end(); ++field) {
    if (field->second->type() == COMPOSITE_VECTOR_FIELD &&
        field->second->io_checkpoint()) {
      bool read_complete = field->second->ReadCheckpoint(checkpoint);
      if (read_complete) field->second->set_initialized();
    }
  }
  
  checkpoint.close_h5file();
  return dt;
};

// Non-member function for checkpointing.
double ReadCheckpointInitialTime(const Comm_ptr_type& comm,
                                 const std::string& filename) {
  HDF5_MPI checkpoint(comm, filename);

  // load the attributes
  double time(0.);
  checkpoint.open_h5file();
  checkpoint.readAttrReal(time, "time");
  checkpoint.close_h5file();
  return time;
};

// Non-member function for checkpointing position.
int ReadCheckpointPosition(const Comm_ptr_type& comm,
                           const std::string& filename) {
  HDF5_MPI checkpoint(comm, filename);

  // load the attributes
  int pos = 0;
  checkpoint.open_h5file();
  checkpoint.readAttrInt(pos, "position");
  checkpoint.close_h5file();
  return pos;
};

// Non-member function for checkpointing observations.
void ReadCheckpointObservations(const Comm_ptr_type& comm,
                                const std::string& filename,
                                Amanzi::ObservationData& obs_data) {
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
  void DeformCheckpointMesh(State& S, Key domain) {
    if (S.HasField(Keys::getKey(domain,"vertex_coordinate"))) { 
      // only deform mesh if vertex coordinate field exists
    AmanziMesh::Mesh * write_access_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S.GetMesh(domain));
    // get vertex coordinates state
    Teuchos::RCP<const CompositeVector> vc = S.GetFieldData(Keys::getKey(domain,"vertex_coordinate"));
    vc->ScatterMasterToGhosted("node");
    const Epetra_MultiVector& vc_n = *vc->ViewComponent("node",true);

    int dim = write_access_mesh_->space_dimension();
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
    if (boost::starts_with(domain, "column"))
      write_access_mesh_->deform( nodeids, new_pos, false, &final_pos); // deforms the mesh
    else
      write_access_mesh_->deform( nodeids, new_pos, true, &final_pos); // deforms the mesh

  }
}


// Non-member function for statistics
void WriteStateStatistics(const State& S, const VerboseObject& vo)
{
  if (vo.os_OK(Teuchos::VERB_HIGH)) {
    Teuchos::OSTab tab = vo.getOSTab();
    *vo.os() << "\nField                                    Min/Max/Avg" << std::endl;

    for (auto f_it = S.field_begin(); f_it != S.field_end(); ++f_it) {
      std::string name(f_it->first);

      if (f_it->second->type() == COMPOSITE_VECTOR_FIELD) {
        std::map<std::string, double> vmin, vmax, vavg;
        f_it->second->GetFieldData()->MinValue(vmin);
        f_it->second->GetFieldData()->MaxValue(vmax);
        f_it->second->GetFieldData()->MeanValue(vavg);

        for (auto c_it = vmin.begin(); c_it != vmin.end(); ++c_it) {
          std::string namedot(name), name_comp(c_it->first);
          if (vmin.size() != 1) namedot.append("." + name_comp);
          namedot.resize(40, '.');
          *vo.os() << namedot << " " << c_it->second << " / " 
                    << vmax[name_comp] << " / " << vavg[name_comp] << std::endl;
        }
      } else if (f_it->second->type() == CONSTANT_SCALAR) {
        double vmin = *f_it->second->GetScalarData();
        name.resize(40, '.');
        *vo.os() << name << " " << vmin << std::endl;
      }
    }
  }
};


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
      Key lifted_key = Keys::getKey(std::get<0>(split)+"_*", std::get<2>(split));
      if (FEList().isParameter(lifted_key)) {
        return FEList().sublist(lifted_key);
      }
    }
  }
  // return an empty new lsit
  return FEList().sublist(key);
}


} // namespace Amanzi
