/*
  State

  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! State, a container for data.
/*!

State  is a  simple data-manager,  allowing PKs  to require,  read, and  write
various fields.

- Acts as a factory for data through the various require methods.
- Provides some data protection by providing both const and non-const
  data pointers to PKs.
- Provides some initialization capability -- this is where all
  independent variables can be initialized (as independent variables
  are owned by state, not by any PK).


.. _state-spec:
.. admonition:: state-spec

   * `"evaluators`" ``[evaluator-typedinline-spec-list]`` A list of evaluators.
      Note this will eventually be an [evaluator-typedinline-spec-list] but the
      evaluators themselves do not include the type info.

   * `"initial conditions`" ``[list]`` A list of constant-in-time variables :
       `"initial conditions`" is a terrible name and will go away in the next
       iteration of state.

.. _evaluator-typedinline-spec:
.. admonition:: evaluator-typedinline-spec

   * `"evaluator type`" ``[string]`` Type of the evaluator Included for
      convenience in defining data that is not in the dependency graph,
      constants are things (like gravity, or atmospheric pressure) which are
      stored in state but never change.  Typically they're limited to scalars
      and dense, local vectors.

.. _constants-scalar-spec:
.. admonition:: constants-scalar-spec

   * `"value`" ``[double]`` Value of a scalar constant

.. _constants-vector-spec:
.. admonition:: constants-vector-spec

   * `"value`" ``[Array(double)]`` Value of a dense, local vector.

Example:

.. code-block:: xml

    <ParameterList name="state">
      <ParameterList name="evaluators">
        <ParameterList name="pressure">
          <Parameter name="evaluator type" type="string" value="primary variable" />
        </ParameterList>
      </ParameterList>

      <ParameterList name="initial conditions">
        <ParameterList name="atmospheric pressure">
          <Parameter name="value" type="double" value="101325.0" />
        </ParameterList>
        <ParameterList name="gravity">
          <Parameter name="value" type="Array(double)" value="{0.0,0.0,-9.80665}" />
        </ParameterList>
      </ParameterList>
    </ParameterList>

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
#include "Checkpoint.hh"
#include "ObservationData.hh"
#include "RecordSet.hh"
#include "StateDefs.hh"
#include "Tag.hh"
#include "Visualization.hh"

namespace Amanzi {

class Evaluator;

enum StateConstructMode {
  STATE_CONSTRUCT_MODE_COPY_POINTERS,
  STATE_CONSTRUCT_MODE_COPY_DATA
};

enum StatePosition {
  TIME_PERIOD_START,
  TIME_PERIOD_INSIDE,
  TIME_PERIOD_END
};

class State {
 private:
  typedef std::map<Key, std::pair<Teuchos::RCP<AmanziMesh::Mesh>, bool> > MeshMap;
  typedef std::map<Key, Teuchos::RCP<AmanziMesh::DomainSet>> DomainSetMap;
  typedef std::map<Key, Teuchos::RCP<Functions::MeshPartition> > MeshPartitionMap;

  using RecordSetMap = std::unordered_map<Key, std::unique_ptr<RecordSet> >;
  using EvaluatorMap = std::unordered_map<Key, std::unordered_map<Tag, Teuchos::RCP<Evaluator> > >;

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
  void InitializeFields();
  void InitializeEvaluators();
  void InitializeFieldCopies();
  bool CheckAllFieldsInitialized();

  // Using another state for initialization (should we use tags? for non-maching states?)
  void Initialize(Teuchos::RCP<State> S);

  // -----------------------------------------------------------------------------
  // State handles mesh management.
  // -----------------------------------------------------------------------------
  // Meshes are "registered" with state. 
  // -- Register a mesh under the default key, "domain".
  void RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                          bool deformable = false);

  // -- Register a mesh under a generic key.
  void RegisterMesh(const Key& key, const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                    bool deformable = false);

  // Ensure a mesh exists.
  bool HasMesh(const Key& key) const { return GetMesh_(key) != Teuchos::null; }
  bool IsDeformableMesh(const Key& key) const;

  // Mesh accessor.
  Teuchos::RCP<const AmanziMesh::Mesh> GetMesh(const Key& key = Key("domain")) const;
  Teuchos::RCP<AmanziMesh::Mesh> GetDeformableMesh(Key key = Key("domain"));

  // Iterate over meshes.
  typedef MeshMap::const_iterator mesh_iterator;
  mesh_iterator mesh_begin() const { return meshes_.begin(); }
  mesh_iterator mesh_end() const { return meshes_.end(); }
  MeshMap::size_type mesh_count() { return meshes_.size(); }

  // DomainSets are collections of meshes, indexed via NAME_GID and referenced
  // to a parent mesh and sets.
  void RegisterDomainSet(const Key& name,
                         const Teuchos::RCP<AmanziMesh::DomainSet> set);
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
  //  F is a factory, which must provide a method Create() that makes a T
  //    (optional)
  template <typename T, typename F>
  F& Require(const Key& fieldname, const Tag& tag, const Key& owner = "") {
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    data_.at(fieldname)->RequireRecord(tag, owner);
    return data_.at(fieldname)->SetType<T, F>();
  }

  template <typename T>
  void Require(const Key& fieldname, const Tag& tag, const Key& owner = "") {
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    data_.at(fieldname)->RequireRecord(tag, owner);
    data_.at(fieldname)->SetType<T>();
  }

  template <typename T, typename F>
  F& Require(const Key& fieldname) {
    return Require<T, F>(fieldname, Tags::DEFAULT, "");
  }

  template <typename T>
  void Require(const Key& fieldname) {
    Require<T>(fieldname, Tags::DEFAULT, "");
  }

  template <typename T, typename F>
  F& Require(const Key& fieldname, const Tag& tag, const Key& owner,
             const std::vector<std::string>& subfield_names) {
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    data_.at(fieldname)->RequireRecord(tag, owner);
    data_.at(fieldname)->GetRecord(tag).set_subfieldnames(subfield_names);
    return data_.at(fieldname)->SetType<T, F>();
  }

  // Ensure a record exists.
  bool HasData(const Key& key, const Tag& tag = Tags::DEFAULT) const {
    if (Keys::hasKey(data_, key)) return data_.at(key)->HasRecord(tag);
    return false;
  }

  // Record accessors.
  const Record& GetRecord(const Key& fieldname, const Tag& tag = Tags::DEFAULT) const {
    return data_.at(fieldname)->GetRecord(tag);
  }
  Record& GetRecordW(const Key& fieldname, const Key& owner) {
    auto& r = data_.at(fieldname)->GetRecord(Tags::DEFAULT);
    r.AssertOwnerOrDie(owner);
    return r;
  }
  Record& GetRecordW(const Key& fieldname, const Tag& tag, const Key& owner) {
    auto& r = data_.at(fieldname)->GetRecord(tag);
    r.AssertOwnerOrDie(owner);
    return r;
  }

  // RecordSet accessors.
  RecordSet& GetRecordSetW(const Key& fieldname) { return *data_.at(fieldname); }

  // Iterate over Records.
  typedef RecordSetMap::const_iterator data_iterator;
  data_iterator data_begin() const { return data_.begin(); }
  data_iterator data_end() const { return data_.end(); }
  RecordSetMap::size_type data_count() { return data_.size(); }

  // Require derivatives
  template <typename T, typename F>
  F& RequireDerivative(const Key& key, const Tag& tag,
                       const Key& wrt_key, const Tag& wrt_tag, const Key& owner = "") {
    auto keytag = Keys::getKeyTag(key, tag.get());
    if (!Keys::hasKey(derivs_, keytag)) {
      derivs_.emplace(keytag, std::make_unique<RecordSet>(keytag));
    }
    Tag dertag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
    derivs_.at(keytag)->RequireRecord(dertag, owner);
    return derivs_.at(keytag)->SetType<T, F>();
  }

  template <typename T>
  void RequireDerivative(const Key& key, const Tag& tag, const Key& wrt_key,
                         const Tag& wrt_tag, const Key& owner = "") {
    auto keytag = Keys::getKeyTag(key, tag.get());
    if (!Keys::hasKey(derivs_, keytag)) {
      derivs_.emplace(keytag, std::make_unique<RecordSet>(keytag));
    }
    Tag dertag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
    derivs_.at(keytag)->RequireRecord(dertag, owner);
    derivs_.at(keytag)->SetType<T>();
  }

  template <typename T, typename F>
  F& RequireDerivative(const Key& key, const Key& wrt_key, const Tag& wrt_tag) {
    return RequireDerivative<T, F>(key, Tags::DEFAULT, wrt_key, wrt_tag, "");
  }

  template <typename T>
  void RequireDerivative(const Key& key, const Key& wrt_key, const Tag& wrt_tag) {
    RequireDerivative<T>(key, Tags::DEFAULT, wrt_key, wrt_tag, "");
  }

  bool HasDerivative(const Key& key, const Tag& tag, const Key& wrt_key, const Tag& wrt_tag) const {
    auto keytag = Keys::getKeyTag(key, tag.get());
    if (Keys::hasKey(derivs_, keytag)) {
      Tag der_tag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
      return derivs_.at(keytag)->HasRecord(der_tag);
    }
    return false;
  }

  bool HasDerivative(const Key& key, const Key& wrt_key) const {
    return HasDerivative(key, Tags::DEFAULT, wrt_key, Tags::DEFAULT);
  }

  // ignoring record access for now, this could be added to, e.g. vis
  // derivatives.
  template <typename T>
  const T& GetDerivative(const Key& key, const Tag& tag, const Key& wrt_key,
                         const Tag& wrt_tag) const {
    Tag der_tag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
    return derivs_.at(Keys::getKeyTag(key, tag.get()))->Get<T>(der_tag);
  }

  template <typename T>
  const T& GetDerivative(const Key& key, const Key& wrt_key) const {
    return GetDerivative<T>(key, Tags::DEFAULT, wrt_key, Tags::DEFAULT);
  }

  template <typename T>
  T& GetDerivativeW(const Key& key, const Tag& tag, const Key& wrt_key,
                    const Tag& wrt_tag, const Key& owner) {
    Tag der_tag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
    return derivs_.at(Keys::getKeyTag(key, tag.get()))->GetW<T>(der_tag, owner);
  }

  template <typename T>
  T& GetDerivativeW(const Key& key, const Key& wrt_key, const Key& owner) {
    return GetDerivativeW<T>(key, Tags::DEFAULT, wrt_key, Tags::DEFAULT, owner);
  }

  template <typename T>
  Teuchos::RCP<const T> GetDerivativePtr(const Key& key, const Tag& tag,
                                         const Key& wrt_key, const Tag& wrt_tag) const {
    Tag der_tag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
    return derivs_.at(Keys::getKeyTag(key, tag.get()))->GetPtr<T>(der_tag);
  }

  template <typename T>
  Teuchos::RCP<T> GetDerivativePtrW(const Key& key, const Tag& tag,
                                    const Key& wrt_key, const Tag& wrt_tag,
                                    const Key& owner) {
    Tag der_tag = make_tag(Keys::getKeyTag(wrt_key, wrt_tag.get()));
    return derivs_.at(Keys::getKeyTag(key, tag.get()))->GetPtrW<T>(der_tag, owner);
  }

  // set operations
  bool HasDerivativeSet(const Key& key, const Tag& tag) const {
    return Keys::hasKey(derivs_, Keys::getKeyTag(key, tag.get()));
  }
  RecordSet& GetDerivativeSet(const Key& key, const Tag& tag) {
    return *derivs_.at(Keys::getKeyTag(key, tag.get()));
  }
  
  // Access to data
  // -- const
  template <typename T>
  const T& Get(const Key& fieldname) const {
    return data_.at(fieldname)->Get<T>();
  }
  template <typename T>
  const T& Get(const Key& fieldname, const Tag& tag) const {
    return data_.at(fieldname)->Get<T>(tag);
  }

  // -- non-const
  template <typename T>
  T& GetW(const Key& fieldname, const Key& owner) {
    return data_.at(fieldname)->GetW<T>(owner);
  }
  template <typename T>
  T& GetW(const Key& fieldname, const Tag& tag, const Key& owner) {
    return data_.at(fieldname)->GetW<T>(tag, owner);
  }

  template <typename T>
  Teuchos::RCP<const T> GetPtr(const Key& fieldname, const Tag& tag = Tags::DEFAULT) const {
    return data_.at(fieldname)->GetPtr<T>(tag);
  }
  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& fieldname, const Tag& tag, const Key& owner) {
    return data_.at(fieldname)->GetPtrW<T>(tag, owner);
  }
  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& fieldname, const Key& owner) {
    return data_.at(fieldname)->GetPtrW<T>(Tags::DEFAULT, owner);
  }

  template <typename T>
  void Set(const Key& fieldname, const Key& owner, const T& data) {
    return data_.at(fieldname)->Set(Tags::DEFAULT, owner, data);
  }
  template <typename T>
  void Set(const Key& fieldname, const Tag& tag, const Key& owner, const T& data) {
    return data_.at(fieldname)->Set(tag, owner, data);
  }

  // -- copy between tags assume that operator = is supported
  template <typename T>
  void Copy(const Key& fieldname, const Tag& tag, const Tag& copy) {
    Key owner = GetRecord(fieldname, copy).owner();
    GetW<T>(fieldname, copy, owner) = Get<T>(fieldname, tag);
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
  Teuchos::ParameterList& GetEvaluatorList(const Key& key);
  
  // -- allows PKs to add to this list to initial conditions
  Teuchos::ParameterList& ICList() { return state_plist_.sublist("initial conditions"); }

  // Evaluator interface
  Evaluator& RequireEvaluator(const Key& key, const Tag& tag = Tags::DEFAULT);

  // -- get/set
  Evaluator& GetEvaluator(const Key& key, const Tag& tag = Tags::DEFAULT);
  const Evaluator& GetEvaluator(const Key& key, const Tag& tag = Tags::DEFAULT) const;
  Teuchos::RCP<Evaluator> GetEvaluatorPtr(const Key& key, const Tag& tag = Tags::DEFAULT);

  void SetEvaluator(const Key& key, const Tag& tag, const Teuchos::RCP<Evaluator>& evaluator);
  void SetEvaluator(const Key& key, const Teuchos::RCP<Evaluator>& evaluator) {
    SetEvaluator(key, Tags::DEFAULT, evaluator);
  }

  bool HasEvaluator(const Key& key, const Tag& tag = Tags::DEFAULT);

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
  Teuchos::ParameterList GetModelParameters(std::string modelname);

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
  double time(const Tag& tag = Tags::DEFAULT) const { return Get<double>("time", tag); }
  void set_time(const Tag& tag, double value) { Set("time", tag, "time", value); }
  void set_time(double value) { Set("time", Tags::DEFAULT, "time", value); }

  void advance_time(const Tag& tag, double dt) {
    Set("time", tag, "time", Get<double>("time", tag) + dt);
  }
  void advance_time(double dt) {
    Set("time", Tags::DEFAULT, "time", Get<double>("time") + dt);
  }

  double final_time() const { return final_time_; }
  void set_final_time(double new_time) { final_time_ = new_time; }
  double intermediate_time() const { return intermediate_time_; }
  void set_intermediate_time(double new_time) { intermediate_time_ = new_time; }

  double last_time() const { return last_time_; }
  void set_last_time( double last_time) { last_time_ = last_time; }
  double initial_time() const { return initial_time_; }
  void set_initial_time( double initial_time) { initial_time_ = initial_time; }

  // Cycle accessor and mutators.
  int cycle() const { return Get<int>("cycle"); }
  void set_cycle(int cycle) { Set("cycle", Tags::DEFAULT, "cycle", cycle); }
  void advance_cycle(int dcycle = 1) {
    Set("cycle", Tags::DEFAULT, "cycle", Get<int>("cycle") + dcycle);
  }

  // Position accessor and mutators.
  int position() const { return Get<int>("position"); }
  void set_position(int pos) { Set("position", Tags::DEFAULT, "position", pos); }

  // Utility for setting vis flags using blacklist and whitelist
  void InitializeIOFlags();

 private:
  // Accessors that return null if the Key does not exist.
  Teuchos::RCP<AmanziMesh::Mesh> GetMesh_(const Key& key) const;
  Teuchos::RCP<const Functions::MeshPartition> GetMeshPartition_(Key);

 private:
  Teuchos::RCP<VerboseObject> vo_;

  // Containers
  MeshMap meshes_;
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

}  // namespace Amanzi

#endif
