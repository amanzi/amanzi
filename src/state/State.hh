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

#ifndef STATE_STATE_HH_
#define STATE_STATE_HH_

// deprecation of old interface
#ifndef STATE_OLD_INTERFACE_DEPRECATED
#  define STATE_OLD_INTERFACE_DEPRECATED 0
#endif

#include <string>
#include <tuple>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "UniqueHelpers.hh"

#if (!STATE_OLD_INTERFACE_DEPRECATED)
#  include "CompositeVector.hh"
#  include "CompositeVectorSpace.hh"
#endif

#include "data/RecordSet.hh"

namespace Amanzi {

class Evaluator;

enum StateConstructMode {
  STATE_CONSTRUCT_MODE_COPY_POINTERS,
  STATE_CONSTRUCT_MODE_COPY_DATA
};
enum StatePosition { TIME_PERIOD_START, TIME_PERIOD_INSIDE, TIME_PERIOD_END };

class State {
 private:
  using MeshMap =
    std::unordered_map<Key, std::pair<Teuchos::RCP<AmanziMesh::Mesh>, bool>>;
  using RecordSetMap = std::unordered_map<Key, std::unique_ptr<RecordSet>>;
  using EvaluatorMap =
    std::unordered_map<Key, std::unordered_map<Key, Teuchos::RCP<Evaluator>>>;
  // using MeshPartitionMap =
  //   std::unordered_map<Key, Teuchos::RCP<Functions::MeshPartition>>;

 public:
  // Default constructor.
  State();

  // Usual constructor.
  explicit State(const Teuchos::RCP<Teuchos::ParameterList>& state_plist);
  explicit State(const Teuchos::ParameterList& state_plist)
    : State(Teuchos::parameterList(state_plist))
  {}

  // // Copy constructor, by default copies memory not pointers.
  // State(const State& other, StateConstructMode
  // mode=STATE_CONSTRUCT_MODE_COPY_DATA);

  // no copy/move constructors, operator=
  State(const State& other) = delete;
  State& operator=(const State&) = delete;
  State(const State&& other) = delete;
  State& operator=(const State&&) = delete;

  // Create data structures, finalizing the structure of the state.
  void Setup();

  // Sub-step in Setup
  void AliasEvaluators();

  // Print state info for debugging
  void Print() const;
  
  void Initialize();

  // -----------------------------------------------------------------------------
  // State handles mesh management.
  // -----------------------------------------------------------------------------
  // Meshes are "registered" with state.  Creation of meshes is NOT handled by
  // state.
  //
  // Register a mesh under the default key, "domain".
  void RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                          bool deformable = false);

  // Register a mesh under a generic key.
  void RegisterMesh(Key key, Teuchos::RCP<AmanziMesh::Mesh> mesh,
                    bool deformable = false);

  // Alias a mesh to an existing mesh
  void AliasMesh(Key target, Key alias);

  // Remove a mesh.
  void RemoveMesh(const Key& key);

  // Ensure a mesh exists.
  bool HasMesh(Key key) const;
  bool IsDeformableMesh(Key key) const;

  // Mesh accessor.
  Teuchos::RCP<const AmanziMesh::Mesh> GetMesh(Key key="domain") const;
  Teuchos::RCP<AmanziMesh::Mesh> GetDeformableMesh(Key key="domain");

  // Iterate over meshes.
  using mesh_iterator = MeshMap::const_iterator;
  mesh_iterator mesh_begin() const { return meshes_.begin(); }
  mesh_iterator mesh_end() const { return meshes_.end(); }
  MeshMap::size_type mesh_count() { return meshes_.size(); }

  // // Some models, typically only defined on cells, are defined by the region.
  // // MeshPartitions are a non-overlapping set of cell regions whose union
  // // covers the mesh.
  // //
  // Teuchos::RCP<const Functions::MeshPartition>
  // GetMeshPartition(const Key& key) const
  // {
  //   return mesh_partitions_.at(key);
  // }

  // Some models are repeated over and over again, especially subgrid models.
  // A flyweight pattern for input parameter lists for such models is enabled.
  // What exactly is done here is open for discussion!
  void RegisterDomainSet(const Key& key) { domain_sets_.insert(key); }

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
  // object -- that object, which is typically either a PK or a
  // Evaluator, may write the solution, and therefore receives non-const
  // pointers to data.  Data may be used by anyone, but non-owning objects
  // receive const-only pointers to data.

  // Require data from State.
  //
  // Takes up to two template arguments:
  //  T is the data type required
  //  F is a factory, which must provide a method Create() that makes a T
  //    (optional)
  template <typename T, typename F>
  F& Require(const Key& fieldname, const Key& tag, const Key& owner = "")
  {
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    RecordSet& r = GetRecordSet(fieldname);
    r.RequireRecord(tag, owner);
    return r.SetType<T, F>();
  }

  template <typename T>
  void Require(const Key& fieldname, const Key& tag, const Key& owner = "")
  {
    if (!Keys::hasKey(data_, fieldname)) {
      data_.emplace(fieldname, std::make_unique<RecordSet>(fieldname));
    }
    RecordSet& r = GetRecordSet(fieldname);
    r.RequireRecord(tag, owner);
    r.SetType<T>();
  }

  template <typename T, typename F>
  F& Require(const Key& fieldname)
  {
    return Require<T, F>(fieldname, "", "");
  }

  template <typename T>
  void Require(const Key& fieldname)
  {
    Require<T>(fieldname, "", "");
  }

  // Ensure a record exists.
  bool HasData(const Key& fieldname, const Key& tag="") const
  {
    if (Keys::hasKey(data_, fieldname)) {
      return GetRecordSet(fieldname).HasRecord(tag);
    }
    return false;
  }

  // Record accessor.
  Record& GetRecordW(const Key& fieldname, const Key& owner)
  {
    Record& r = GetRecordSet(fieldname).GetRecord("");
    r.AssertOwnerOrDie(owner);
    return r;
  }
  Record& GetRecordW(const Key& fieldname, const Key& tag, const Key& owner)
  {
    Record& r = GetRecordSet(fieldname).GetRecord(tag);
    r.AssertOwnerOrDie(owner);
    return r;
  }
  const Record& GetRecord(const Key& fieldname, const Key& tag = "") const
  {
    return GetRecordSet(fieldname).GetRecord(tag);
  }

  // Record accessor.
  RecordSet& GetRecordSet(const Key& fieldname);
  const RecordSet& GetRecordSet(const Key& fieldname) const;
  
  // Iterate over Records.
  typedef RecordSetMap::const_iterator data_iterator;
  data_iterator data_begin() const { return data_.begin(); }
  data_iterator data_end() const { return data_.end(); }
  RecordSetMap::size_type data_count() { return data_.size(); }

  // Access to data
  template <typename T>
  const T& Get(const Key& fieldname) const
  {
    return GetRecordSet(fieldname).Get<T>();
  }
  template <typename T>
  const T& Get(const Key& fieldname, const Key& tag) const
  {
    return GetRecordSet(fieldname).Get<T>(tag);
  }
  template <typename T>
  T& GetW(const Key& fieldname, const Key& owner)
  {
    return GetRecordSet(fieldname).GetW<T>(owner);
  }
  template <typename T>
  T& GetW(const Key& fieldname, const Key& tag, const Key& owner)
  {
    return GetRecordSet(fieldname).GetW<T>(tag, owner);
  }
  template <typename T>
  Teuchos::RCP<const T> GetPtr(const Key& fieldname, const Key& tag = "") const
  {
    return GetRecordSet(fieldname).GetPtr<T>(tag);
  }
  template <typename T>
  Teuchos::RCP<T> GetPtrW(const Key& fieldname, const Key& owner)
  {
    return GetRecordSet(fieldname).GetPtrW<T>("", owner);
  }
  template <typename T>
  Teuchos::RCP<T>
  GetPtrW(const Key& fieldname, const Key& tag, const Key& owner)
  {
    return GetRecordSet(fieldname).GetPtrW<T>(tag, owner);
  }

  template <typename T>
  void Set(const Key& fieldname, const Key& owner, const T& data)
  {
    GetRecordSet(fieldname).Set("", owner, data);
  }
  template <typename T>
  void
  Set(const Key& fieldname, const Key& tag, const Key& owner, const T& data)
  {
    GetRecordSet(fieldname).Set(tag, owner, data);
  }
  template <typename T>
  void
  SetPtr(const Key& fieldname, const Key& owner, const Teuchos::RCP<T>& data)
  {
    GetRecordSet(fieldname).SetPtr("", owner, data);
  }
  template <typename T>
  void SetPtr(const Key& fieldname, const Key& tag, const Key& owner,
              const Teuchos::RCP<T>& data)
  {
    GetRecordSet(fieldname).SetPtr(tag, owner, data);
  }


  template <typename T, typename F>
  F& RequireDerivative(const Key& key, const Key& tag, const Key& wrt_key,
                       const Key& wrt_tag, const Key& owner = "")
  {
    auto keytag = Keys::getKeyTag(key, tag);
    if (!Keys::hasKey(derivs_, keytag)) {
      derivs_.emplace(keytag, std::make_unique<RecordSet>(keytag));
    }
    derivs_.at(keytag)->RequireRecord(Keys::getKeyTag(wrt_key, wrt_tag), owner);
    return derivs_.at(keytag)->SetType<T, F>();
  }

  template <typename T>
  void RequireDerivative(const Key& key, const Key& tag, const Key& wrt_key,
                         const Key& wrt_tag, const Key& owner = "")
  {
    auto keytag = Keys::getKeyTag(key, tag);
    if (!Keys::hasKey(derivs_, keytag)) {
      derivs_.emplace(keytag, std::make_unique<RecordSet>(keytag));
    }
    derivs_.at(keytag)->RequireRecord(Keys::getKeyTag(wrt_key, wrt_tag), owner);
    derivs_.at(keytag)->SetType<T>();
  }

  template <typename T, typename F>
  F& RequireDerivative(const Key& key, const Key& wrt_key, const Key& wrt_tag)
  {
    return RequireDerivative<T, F>(key, "", wrt_key, wrt_tag, "");
  }

  template <typename T>
  void RequireDerivative(const Key& key, const Key& wrt_key, const Key& wrt_tag)
  {
    RequireDerivative<T>(key, "", wrt_key, wrt_tag, "");
  }

  bool HasDerivative(const Key& key, const Key& tag, const Key& wrt_key,
                     const Key& wrt_tag) const
  {
    auto keytag = Keys::getKeyTag(key, tag);
    if (Keys::hasKey(derivs_, keytag)) {
      return derivs_.at(keytag)->HasRecord(Keys::getKeyTag(wrt_key, wrt_tag));
    }
    return false;
  }

  // ignoring record access for now, this could be added to, e.g. vis
  // derivatives.

  template <typename T>
  const T& GetDerivative(const Key& key, const Key& tag, const Key& wrt_key,
                         const Key& wrt_tag) const
  {
    return GetDerivativeSet(key, tag).Get<T>(Keys::getKeyTag(wrt_key, wrt_tag));
  }

  template <typename T>
  T& GetDerivativeW(const Key& key, const Key& tag, const Key& wrt_key,
                    const Key& wrt_tag, const Key& owner)
  {
    return GetDerivativeSet(key, tag).GetW<T>(Keys::getKeyTag(wrt_key, wrt_tag), owner);
  }

  template <typename T>
  Teuchos::RCP<const T>
  GetDerivativePtr(const Key& key, const Key& tag, const Key& wrt_key,
                   const Key& wrt_tag) const
  {
    return GetDerivativeSet(key, tag).GetPtr<T>(Keys::getKeyTag(wrt_key, wrt_tag));
  }

  template <typename T>
  Teuchos::RCP<T>
  GetDerivativePtrW(const Key& key, const Key& tag, const Key& wrt_key,
                    const Key& wrt_tag, const Key& owner)
  {
    return GetDerivativeSet(key, tag).GetPtrW<T>(Keys::getKeyTag(wrt_key, wrt_tag), owner);
  }


  bool HasDerivativeSet(const Key& key, const Key& tag) const
  {
    return Keys::hasKey(derivs_, Keys::getKeyTag(key, tag));
  }
  RecordSet& GetDerivativeSet(const Key& key, const Key& tag);
  const RecordSet& GetDerivativeSet(const Key& key, const Key& tag) const;

  // A few special parameters with special methods
  double time(const Key& tag = "") const { return Get<double>("time", tag); }
  void set_time(const Key& tag, double value)
  {
    Set("time", tag, "time", value);
  }
  void set_time(double value) { Set("time", "", "time", value); }
  void advance_time(const Key& tag, double dt)
  {
    Set("time", tag, "time", Get<double>("time", tag) + dt);
  }
  void advance_time(double dt)
  {
    Set("time", "", "time", Get<double>("time") + dt);
  }

  int cycle(const Key& tag = "") const { return Get<int>("cycle", tag); }
  void set_cycle(const Key& tag, int value)
  {
    Set("cycle", tag, "cycle", value);
  }
  void set_cycle(int value) { Set("cycle", "", "cycle", value); }
  void advance_cycle(const Key& tag, int dc = 1)
  {
    Set("cycle", tag, "cycle", Get<int>("cycle", tag) + dc);
  }
  void advance_cycle(int dc = 1)
  {
    Set("cycle", "", "cycle", Get<int>("cycle") + dc);
  }

  // -----------------------------------------------------------------------------
  // State handles data evaluation.
  // -----------------------------------------------------------------------------
  // To manage lazy yet sufficient updating of models and derivatives of
  // models, we use a graph-based view of data and data dependencies, much
  // like the Phalanx approach.  A directed acyclic graph of dependencies are
  // managed in State, where each node is a Evaluator.
  //
  // Access to the FEList -- this allows PKs to add to this list for custom
  // evaluators.
  Teuchos::ParameterList& FEList()
  {
    return state_plist_->sublist("evaluators");
  }

  // Require Evaluators.
  Evaluator& RequireEvaluator(const Key&, const Key& tag = "");

  // Ensure a Evaluator exists.
  bool HasEvaluator(const Key&, const Key& tag = "");

  // Evaluator accessor.
  Evaluator& GetEvaluator(const Key&, const Key& tag = "");
  const Evaluator& GetEvaluator(const Key&, const Key& tag = "") const;
  Teuchos::RCP<Evaluator> GetEvaluatorPtr(const Key&, const Key& tag = "");

  // Evaluator mutator.
  void SetEvaluator(const Key& key, const Teuchos::RCP<Evaluator>& evaluator)
  {
    SetEvaluator(key, "", evaluator);
  }
  void SetEvaluator(const Key& key, const Key& tag,
                    const Teuchos::RCP<Evaluator>& evaluator);

  // Iterate over evaluators.
  typedef EvaluatorMap::const_iterator evaluator_iterator;
  evaluator_iterator evaluator_begin() const { return evaluators_.begin(); }
  evaluator_iterator evaluator_end() const { return evaluators_.end(); }
  EvaluatorMap::size_type evaluator_count() { return evaluators_.size(); }

  // Write evaluators to file for drawing dependency graph.
  void WriteDependencyGraph() const;
  // void WriteStatistics(Teuchos::RCP<VerboseObject>& vo) const;

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

 private:
  // Iterate over Records.
  typedef RecordSetMap::iterator data_iterator_;
  data_iterator_ data_begin_() { return data_.begin(); }
  data_iterator_ data_end_() { return data_.end(); }

  // read checkpoint is friend to have access to non-const iterator
  friend void
  ReadCheckpoint(const Comm_ptr_type& comm, State&, const std::string&);

  // Containers
  MeshMap meshes_;
  RecordSetMap data_;
  RecordSetMap derivs_;
  EvaluatorMap evaluators_;
  // MeshPartitionMap mesh_partitions_;
  std::set<Key> domain_sets_;

  // parameter list
  Teuchos::RCP<Teuchos::ParameterList> state_plist_;
  Teuchos::RCP<VerboseObject> vo_;
};

// -----------------------------------------------------------------------------
// Non-member functions for I/O of a State.
// -----------------------------------------------------------------------------

// Visualization of State.
void
WriteVis(Visualization& vis, State& S, const Key& tag="next");

// Checkpointing State.
void
WriteCheckpoint(Checkpoint& chkp, const State& S, bool final = false);

void
ReadCheckpoint(const Comm_ptr_type& comm, State& S,
               const std::string& filename);

void
DeformCheckpointMesh(State& S);

} // namespace Amanzi

#endif
