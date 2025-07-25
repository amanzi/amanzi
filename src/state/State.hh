/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

State is a simple data manager, allowing PKs to require, read, and write
various fields. It:

- Acts as a factory for data through the various require methods.
- Provides some data protection by providing both const and non-const
  access to data.
- Provides some initialization capability.
- Can be visualized or checkpointed.

.. _state-spec:
.. admonition:: state-spec

   * `"evaluators`" ``[evaluator-typedinline-spec-list]`` A list of evaluators.

   * `"initial conditions`" ``[list]`` A list of constant-in-time data and
     other constants (e.g. gravity and atmospheric pressure).  Note that
     `"initial conditions`" is not a particularly descriptive name here -- PK
     initial conditions are generally not here.

   * `"model parameters`" ``[list]`` A list of shared model parameters that can
     be used across all evaluators.

   * `"debug`" ``[state-debug-spec]`` **optional**


State also houses the dependency graph of all evaluators.

Evaluators are individual terms used to build up a PK or MPCs.  Each
term represents a variable in the equation, and can consist of primary
variables (those that are solved for by a PK solver), independent
variables (those that have no dependent variables but are provided by
the user as data), and secondary variables (those that are functions
of other variables).  Note that all three may be variable in space
and/or time.


.. _evaluator-typedinline-spec:
.. admonition:: evaluator-typedinline-spec

   * `"evaluator type`" ``[string]`` Type of the evaluator


State's debug spec is intended for use by PK and Evaluator developers, or occasionally
users debugging input files.

.. _state-debug-spec:
.. admonition:: state-debug-spec

   * `"data`" ``[Array(string)]`` A list of data.  If a variable is in this
     list, a message is printed every time that data is required.  Setting a
     breakpoint in a debugger on a line of code in
     `State::CheckIsDebugData_()` will allow debugging all locations that
     require this data, enabling debugging of incorrect or inconsistent require
     statements.

   * `"evaluators`" ``[Array(string)]`` A list of evaluators.  If a variable is
     in this list, a message is printed every time an evaluator for this
     variable is reuqired.  Setting a breakpoint in a debugger in
     `State::CheckIsDebugEval_()` will allow debugging of all locations that
     require a given evaluator.

*/

#ifndef STATE_STATE_HH_
#define STATE_STATE_HH_

#include <string>
#include <vector>

// TPLs
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

// Amanzi
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "DomainSet.hh"
#include "Key.hh"
#include "Mesh.hh"
#include "MeshPartition.hh"

// Amanzi::State
#include "ObservationData.hh"
#include "RecordSet.hh"
#include "StateDefs.hh"
#include "Tag.hh"

namespace Amanzi {

class Evaluator;
class EvaluatorAlias;

enum StateConstructMode { STATE_CONSTRUCT_MODE_COPY_POINTERS, STATE_CONSTRUCT_MODE_COPY_DATA };

enum StatePosition { TIME_PERIOD_START, TIME_PERIOD_INSIDE, TIME_PERIOD_END };

class State {
 private:
  typedef std::map<Key, std::pair<Teuchos::RCP<AmanziMesh::Mesh>, bool>> MeshMap;
  typedef std::map<Key, Teuchos::RCP<AmanziMesh::DomainSet>> DomainSetMap;
  typedef std::map<Key, Teuchos::RCP<Functions::MeshPartition>> MeshPartitionMap;

  using RecordSetMap = std::map<Key, std::unique_ptr<RecordSet>>;
  using EvaluatorMap = std::map<Key, std::unordered_map<Tag, Teuchos::RCP<Evaluator>>>;

 public:
  State();
  explicit State(Teuchos::ParameterList& state_plist);

  // Copy constructor, copies memory not pointers.
  // State(const State& other, StateConstructMode mode=STATE_CONSTRUCT_MODE_COPY_DATA);

  // Assignment and copy operators. Note this should be replaced with smart
  // usage of tags
  State(const State& other) = delete;
  State& operator=(const State& other) = delete;
  State(const State&& other) = delete;
  State& operator=(const State&& other) = delete;

  // Set requirements from all evaluators, calling EnsureCompatibility and
  // allocating all memory.
  void Setup();

  // Sub-steps in the initialization process. (Used by Amanzi)
  void Initialize();
  void InitializeFields(const Tag& tag = Tags::DEFAULT);
  void InitializeEvaluators();
  void InitializeFieldCopies(const Tag& ref = Tags::DEFAULT);
  bool CheckAllFieldsInitialized();

  // Using another state for initialization (should we use tags? for
  // non-maching states?)
  void Initialize(const State& S);

  // -----------------------------------------------------------------------------
  // State handles mesh management.
  // -----------------------------------------------------------------------------
  // Meshes are "registered" with state.
  // -- Register a mesh under the default key, "domain".
  void RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh, bool deformable = false);

  // -- Register a mesh under a generic key.
  void RegisterMesh(const Key& key,
                    const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                    bool deformable = false);

  // Alias a mesh to an existing mesh
  void AliasMesh(const Key& target, const Key& alias);
  bool IsAliasedMesh(const Key& key) const;

  // Ensure a mesh exists.
  bool HasMesh(const Key& key) const { return GetMesh_(key) != Teuchos::null; }
  bool IsDeformableMesh(const Key& key) const;

  // Mesh accessor.
  Teuchos::RCP<AmanziMesh::Mesh> GetMesh(const Key& key = Key("domain")) const;
  Teuchos::RCP<AmanziMesh::Mesh> GetDeformableMesh(Key key = Key("domain"));

  // Iterate over meshes.
  typedef MeshMap::const_iterator mesh_iterator;
  mesh_iterator mesh_begin() const { return meshes_.begin(); }
  mesh_iterator mesh_end() const { return meshes_.end(); }
  MeshMap::size_type mesh_count() { return meshes_.size(); }

  // DomainSets are collections of meshes, indexed via NAME_GID and referenced
  // to a parent mesh and sets.
  void RegisterDomainSet(const Key& name, const Teuchos::RCP<AmanziMesh::DomainSet> set);
  bool HasDomainSet(const Key& name) const;
  Teuchos::RCP<const AmanziMesh::DomainSet> GetDomainSet(const Key& name) const;

  // -----------------------------------------------------------------------------
  // State handles data management.
  // -----------------------------------------------------------------------------
  // Data is managed by a Record, which both controls access and provides
  // metadata.
  //
  // State manages the creation and consistency of data.  Data is "required"
  // of the state.  The require methods act as factories and consistency
  // checks for ownership and type specifiers of the fields.
  //
  // State also manages access to data.  Data is "owned" by at most one
  // object -- that object, which is typically either a PK or a Evaluator
  // (always true for secondary), may write the solution, and therefore
  // receives non-const pointers to data.  Data may be used by anyone, but
  // non-owning objects receive const-only pointers to data.
  //
  // Requiring data from State takes up to two template arguments:
  //  T is the data type required
  //  F is a factory, which must provide a method Create() that makes a T, or
  //    can be used in a constructor, e.g. T(F), or is optional if T needs no
  //    metadata to be created, e.g. T == double.
  //
  //
  // This Require call will not compile for factories F that do not have a
  // default constructor (e.g. Epetra_Map).
  //
  // Note, aliases are created for tag == NEXT when:
  //  1. there is no record at NEXT
  //  2. there _is_ a record at another tag including next in the name
  //  3. alias_ok is true
  // This allows e.g. observations to use aliases for variables that would not
  // otherwise exist, because they are only used in subcycled PKs, rather than
  // forcing the creation and update of a true NEXT variable.
  template<typename T, typename F>
  F& Require(const Key& fieldname, const Tag& tag, const Key& owner = "", bool alias_ok = false)
  {
    AMANZI_ASSERT(fieldname != "");
    CheckIsDebugData_(fieldname, tag);
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    GetRecordSetW(fieldname).RequireRecord(tag, owner, alias_ok);
    return GetRecordSetW(fieldname).SetType<T, F>();
  }

  // This Require call is for factories that do not have a default constructor.
  // (e.g. Epetra_Map).
  template<typename T, typename F>
  F& Require(const F& f,
             const Key& fieldname,
             const Tag& tag,
             const Key& owner = "",
             bool alias_ok = true)
  {
    AMANZI_ASSERT(fieldname != "");
    CheckIsDebugData_(fieldname, tag);
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    GetRecordSetW(fieldname).RequireRecord(tag, owner, alias_ok);
    return GetRecordSetW(fieldname).SetType<T, F>(f);
  }

  // This Require call is for default-constructible data that does not need a
  // factory (e.g. plain-old-data like int, double, etc).
  template<typename T>
  void Require(const Key& fieldname, const Tag& tag, const Key& owner = "", bool alias_ok = true)
  {
    AMANZI_ASSERT(fieldname != "");
    CheckIsDebugData_(fieldname, tag);
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    auto& rs = GetRecordSetW(fieldname);
    rs.RequireRecord(tag, owner, alias_ok);
    rs.SetType<T>();
  }

  // This Require call will not compile for factories F that do not have a
  // default constructor (e.g. Epetra_Map).
  template<typename T, typename F>
  F& Require(const Key& fieldname,
             const Tag& tag,
             const Key& owner,
             const std::vector<std::string>& subfield_names,
             bool alias_ok = true)
  {
    AMANZI_ASSERT(fieldname != "");
    CheckIsDebugData_(fieldname, tag);
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    auto& rs = GetRecordSetW(fieldname);
    rs.set_subfieldnames(subfield_names);
    rs.RequireRecord(tag, owner, alias_ok);
    return rs.SetType<T, F>();
  }

  // This Require call is for factories that do not have a default constructor.
  // (e.g. Epetra_Map).
  template<typename T, typename F>
  F& Require(const F& f,
             const Key& fieldname,
             const Tag& tag,
             const Key& owner,
             const std::vector<std::string>& subfield_names,
             bool alias_ok = true)
  {
    AMANZI_ASSERT(fieldname != "");
    CheckIsDebugData_(fieldname, tag);
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    auto& rs = GetRecordSetW(fieldname);
    rs.set_subfieldnames(subfield_names);
    rs.RequireRecord(tag, owner, alias_ok);
    return rs.SetType<T, F>(f);
  }

  //
  // RecordSets are a collection of Records that share the same data factory
  // (and therefore data layout).
  //
  bool HasRecordSet(const Key& fieldname) const;
  const RecordSet& GetRecordSet(const Key& fieldname) const;
  RecordSet& GetRecordSetW(const Key& fieldname);

  // Iterate over RecordSets
  typedef RecordSetMap::const_iterator data_iterator;
  data_iterator data_begin() const { return data_.begin(); }
  data_iterator data_end() const { return data_.end(); }
  RecordSetMap::size_type data_count() { return data_.size(); }

  //
  // Records are an instance of data + metadata, created by the RecordSet's
  // data factory, and distinguished by tags.
  //
#ifndef DISABLE_DEFAULT_TAG
  bool HasRecord(const Key& key, const Tag& tag = Tags::DEFAULT) const;
  const Record& GetRecord(const Key& fieldname, const Tag& tag = Tags::DEFAULT) const;
#else
  bool HasRecord(const Key& key, const Tag& tag) const;
  const Record& GetRecord(const Key& fieldname, const Tag& tag) const;
#endif

  Record& GetRecordW(const Key& fieldname, const Key& owner);
  Record& GetRecordW(const Key& fieldname, const Tag& tag, const Key& owner);

  // Require derivatives
  template<typename T, typename F>
  F& RequireDerivative(const Key& key,
                       const Tag& tag,
                       const Key& wrt_key,
                       const Tag& wrt_tag,
                       const Key& owner = "")
  {
    if (!HasDerivativeSet(key, tag)) {
      auto keytag = Keys::getKey(key, tag);
      derivs_.emplace(keytag, std::make_unique<RecordSet>(keytag));
    }
    Tag dertag(Keys::getKey(wrt_key, wrt_tag));
    auto& deriv_set = GetDerivativeSetW(key, tag);
    deriv_set.RequireRecord(dertag, owner);
    return deriv_set.SetType<T, F>();
  }

  template<typename T>
  void RequireDerivative(const Key& key,
                         const Tag& tag,
                         const Key& wrt_key,
                         const Tag& wrt_tag,
                         const Key& owner = "")
  {
    if (!HasDerivativeSet(key, tag)) {
      auto keytag = Keys::getKey(key, tag);
      derivs_.emplace(keytag, std::make_unique<RecordSet>(keytag));
    }
    Tag dertag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    auto& deriv_set = GetDerivativeSetW(key, tag);
    deriv_set.RequireRecord(dertag, owner);
    return deriv_set.SetType<T>();
  }

  // set operations
  bool HasDerivativeSet(const Key& key, const Tag& tag) const;
  const RecordSet& GetDerivativeSet(const Key& key, const Tag& tag) const;
  RecordSet& GetDerivativeSetW(const Key& key, const Tag& tag);

  bool HasDerivative(const Key& key, const Tag& tag, const Key& wrt_key, const Tag& wrt_tag) const;
  bool HasDerivative(const Key& key, const Key& wrt_key) const;

  // ignoring record access for now, this could be added to, e.g. vis
  // derivatives.
  template<typename T>
  const T& GetDerivative(const Key& key,
                         const Tag& tag,
                         const Key& wrt_key,
                         const Tag& wrt_tag) const
  {
    Tag der_tag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    return GetDerivativeSet(key, tag).Get<T>(der_tag);
  }

  template<typename T>
  T& GetDerivativeW(const Key& key,
                    const Tag& tag,
                    const Key& wrt_key,
                    const Tag& wrt_tag,
                    const Key& owner)
  {
    Tag der_tag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    return GetDerivativeSetW(key, tag).GetW<T>(der_tag, owner);
  }

  template<typename T>
  Teuchos::RCP<const T> GetDerivativePtr(const Key& key,
                                         const Tag& tag,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag) const
  {
    Tag der_tag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    return GetDerivativeSet(key, tag).GetPtr<T>(der_tag);
  }

  template<typename T>
  Teuchos::RCP<T> GetDerivativePtrW(const Key& key,
                                    const Tag& tag,
                                    const Key& wrt_key,
                                    const Tag& wrt_tag,
                                    const Key& owner)
  {
    Tag der_tag = make_tag(Keys::getKey(wrt_key, wrt_tag));
    return GetDerivativeSetW(key, tag).GetPtrW<T>(der_tag, owner);
  }

  // Access to data
#ifndef DISABLE_DEFAULT_TAG
  template<typename T>
  bool HasData(const Key& fieldname, const Tag& tag = Tags::DEFAULT) const
  {
    return HasRecord(fieldname, tag) && GetRecord(fieldname, tag).ValidType<T>();
  }

  // -- const
  template<typename T>
  const T& Get(const Key& fieldname, const Tag& tag = Tags::DEFAULT) const
  {
    return GetRecordSet(fieldname).Get<T>(tag);
  }

  // -- non-const
  template<typename T>
  T& GetW(const Key& fieldname, const Key& owner)
  {
    return GetW<T>(fieldname, Tags::DEFAULT, owner);
  }

#else

  template<typename T>
  bool HasData(const Key& fieldname, const Tag& tag) const
  {
    return HasRecord(fieldname, tag) && GetRecord(fieldname, tag).ValidType<T>();
  }

  // -- const
  template<typename T>
  const T& Get(const Key& fieldname, const Tag& tag) const
  {
    return GetRecordSet(fieldname).Get<T>(tag);
  }

#endif

  template<typename T>
  T& GetW(const Key& fieldname, const Tag& tag, const Key& owner)
  {
    return GetRecordSetW(fieldname).GetW<T>(tag, owner);
  }

  template<typename T>
  Teuchos::RCP<const T> GetPtr(const Key& fieldname, const Tag& tag) const
  {
    return GetRecordSet(fieldname).GetPtr<T>(tag);
  }

  template<typename T>
  Teuchos::RCP<T> GetPtrW(const Key& fieldname, const Tag& tag, const Key& owner)
  {
    return GetRecordSetW(fieldname).GetPtrW<T>(tag, owner);
  }

  //
  // Sets by deep copy, not pointer
  //
  template<typename T>
  void Assign(const Key& fieldname, const Tag& tag, const Key& owner, const T& data)
  {
    return GetRecordSetW(fieldname).Assign(tag, owner, data);
  }


  // Sets by pointer
  template<typename T>
  void SetPtr(const Key& fieldname, const Tag& tag, const Key& owner, const Teuchos::RCP<T>& data)
  {
    return GetRecordSetW(fieldname).SetPtr<T>(tag, owner, data);
  }

  // -- assign between tags assume that operator= on the data type is supported
  inline void Assign(const Key& fieldname, const Tag& dest, const Tag& source)
  {
    return GetRecordSetW(fieldname).Assign(dest, source);
  }


  // -----------------------------------------------------------------------------
  // State handles data evaluation.
  // -----------------------------------------------------------------------------
  // To manage lazy yet sufficient updating of models and derivatives of
  // models, we use a graph-based view of data and data dependencies, much
  // like the Phalanx approach.  A directed acyclic graph of dependencies are
  // managed in State, where each node is an Evaluator.
  //
  // -- allows PKs to add to this list to custom evaluators
  Teuchos::ParameterList& FEList() { return state_plist_.sublist("evaluators"); }
  const Teuchos::ParameterList& FEList() const { return state_plist_.sublist("evaluators"); }
  Teuchos::ParameterList& GetEvaluatorList(const Key& key);
  bool HasEvaluatorList(const Key& key) const;

  // -- allows PKs to add to this list to initial conditions
  Teuchos::ParameterList& ICList() { return state_plist_.sublist("initial conditions"); }
  const Teuchos::ParameterList& ICList() const
  {
    return state_plist_.sublist("initial conditions");
  }
  Teuchos::ParameterList& GetICList(const Key& key);
  bool HasICList(const Key& key) const;

  // Evaluator interface
  Evaluator& RequireEvaluator(const Key& key, const Tag& tag, bool alias_ok = true);

#ifndef DISABLE_DEFAULT_TAG
  // -- get/set
  Evaluator& GetEvaluator(const Key& key, const Tag& tag = Tags::DEFAULT);
  const Evaluator& GetEvaluator(const Key& key, const Tag& tag = Tags::DEFAULT) const;
#else
  // -- get/set
  Evaluator& GetEvaluator(const Key& key, const Tag& tag);
  const Evaluator& GetEvaluator(const Key& key, const Tag& tag) const;
#endif

  bool HasEvaluator(const Key& key, const Tag& tag) const;
  void SetEvaluator(const Key& key, const Tag& tag, const Teuchos::RCP<Evaluator>& evaluator);
  Teuchos::RCP<Evaluator> GetEvaluatorPtr(const Key& key, const Tag& tag);

  // -- iterators/counts
  int evaluator_count() { return evaluators_.size(); }

  // Write evaluators to file for drawing dependency graph.
  void WriteDependencyGraph() const;

  // -----------------------------------------------------------------------------
  // State handles model parameters.
  // -----------------------------------------------------------------------------
  // Some model parameters may be common to many PKs, Evaluators, boundary
  // conditions, etc.  Access to the parameters required to make these models
  // is handled through state.  This is used infrequently currently, and
  // should be used and tested more thoroughly.
  //
  // Get a parameter list.
  const Teuchos::ParameterList& GetModelParameters(const std::string& modelname) const;

  // -----------------------------------------------------------------------------
  // State handles MeshPartitions
  // -----------------------------------------------------------------------------
  // Some models, typically only defined on cells, are defined by the region.
  // MeshPartitions are a non-overlapping set of cell regions whose union
  // covers the mesh.
  //
  Teuchos::RCP<const Functions::MeshPartition> GetMeshPartition(Key);

  // -----------------------------------------------------------------------------
  // Time tags and vector copies
  // -----------------------------------------------------------------------------
  // Time accessor and mutators.
  void require_time(const Tag& tag, const Key& owner = "time");
  double get_time(const Tag& tag = Tags::DEFAULT) const { return Get<double>("time", tag); }
  void set_time(const Tag& tag, double value);
  void set_time(double value) { Assign("time", Tags::DEFAULT, "time", value); }

  void advance_time(const Tag& tag, double dt) { set_time(tag, get_time(tag) + dt); }
  void advance_time(double dt) { advance_time(Tags::DEFAULT, dt); }

  // can these go away in favor of time at different tags?
  double final_time() const { return final_time_; }
  void set_final_time(double new_time) { final_time_ = new_time; }
  double intermediate_time() const { return intermediate_time_; }
  void set_intermediate_time(double new_time) { intermediate_time_ = new_time; }

  double last_time() const { return last_time_; }
  void set_last_time(double last_time) { last_time_ = last_time; }
  double initial_time() const { return initial_time_; }
  void set_initial_time(double initial_time) { initial_time_ = initial_time; }

  // Cycle accessor and mutators.
  void require_cycle(const Tag& tag) { Require<int>("cycle", tag, "cycle", false); }
  int get_cycle(Tag tag = Tags::DEFAULT) const { return Get<int>("cycle", tag); }
  void set_cycle(Tag tag, int cycle) { Assign("cycle", tag, "cycle", cycle); }
  void set_cycle(int cycle) { set_cycle(Tags::DEFAULT, cycle); }
  void advance_cycle(Tag tag = Tags::DEFAULT, int dcycle = 1)
  {
    Assign("cycle", tag, "cycle", get_cycle(tag) + dcycle);
  }

  // Position accessor and mutators.
  int get_position() const { return Get<int>("position", Tags::DEFAULT); }
  void set_position(int pos) { Assign("position", Tags::DEFAULT, "position", pos); }

  // Utility for setting vis flags using blacklist and whitelist
  void InitializeIOFlags();

 private:
  // Accessors that return null if the Key does not exist.
  Teuchos::RCP<AmanziMesh::Mesh> GetMesh_(const Key& key) const;
  Teuchos::RCP<const Functions::MeshPartition> GetMeshPartition_(Key);
  void SetEvaluator_(const Key& key, const Tag& tag, const Teuchos::RCP<Evaluator>& evaluator);

  // a hook to allow debuggers to connect
  bool CheckIsDebugEval_(const Key& key, const Tag& tag);
  bool CheckIsDebugData_(const Key& key, const Tag& tag);

 private:
  Teuchos::RCP<VerboseObject> vo_;

  // Containers
  MeshMap meshes_;
  std::map<Key, Key> mesh_aliases_;
  RecordSetMap data_;
  RecordSetMap derivs_;
  EvaluatorMap evaluators_;

  MeshPartitionMap mesh_partitions_;
  DomainSetMap domain_sets_;

  // meta-data
  double final_time_, intermediate_time_, last_time_, initial_time_;

  // parameter list
  Teuchos::ParameterList state_plist_;
};

} // namespace Amanzi

#endif
