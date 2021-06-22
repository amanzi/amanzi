/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
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

   * `"field evaluators`" ``[field-evaluator-typedinline-spec-list]`` A list of evaluators.
      Note this will eventually be an [evaluator-typedinline-spec-list] but the
      evaluators themselves do not include the type info.

   * `"initial conditions`" ``[list]`` A list of constant-in-time variables :
       `"initial conditions`" is a terrible name and will go away in the next
       iteration of state.

.. _field-evaluator-typedinline-spec:
.. admonition:: field-evaluator-typedinline-spec

   * `"field evaluator type`" ``[string]`` Type of the evaluator Included for
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
      <ParameterList name="field evaluators">
        <ParameterList name="pressure">
          <Parameter name="field evaluator type" type="string" value="primary variable field evaluator" />
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
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"
#include "DomainSet.hh"
#include "MeshPartition.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"

#include "Visualization.hh"
#include "Checkpoint.hh"
#include "ObservationData.hh"

#include "Field.hh"
#include "Field_Scalar.hh"
#include "Field_ConstantVector.hh"
#include "Field_CompositeVector.hh"

namespace Amanzi {

class FieldEvaluator;

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

  typedef std::map<Key, std::pair<Teuchos::RCP<AmanziMesh::Mesh>,
                                  bool> > MeshMap;
  typedef std::map<Key, Teuchos::RCP<CompositeVectorSpace> > FieldFactoryMap;
  typedef std::map<Key, Teuchos::RCP<Field> > FieldMap;
  typedef std::map<Key, Teuchos::RCP<FieldEvaluator> > FieldEvaluatorMap;
  typedef std::map<Key, Teuchos::RCP<AmanziMesh::DomainSet>> DomainSetMap;
  typedef std::map<Key, Teuchos::RCP<Functions::MeshPartition> > MeshPartitionMap;

 public:

  // Default constructor.
  State();

  // Usual constructor.
  explicit State(Teuchos::ParameterList& state_plist);

  // Copy constructor, copies memory not pointers.
  State(const State& other, StateConstructMode mode=STATE_CONSTRUCT_MODE_COPY_DATA);

  // Assignment operator, copies memory not pointers.  Note this
  // implementation requires the State being copied has the same structure (in
  // terms of fields, order of fields, etc) as *this.  This really means that
  // it should be a previously-copy-constructed version of the State.  One and
  // only one State should be instantiated and populated -- all other States
  // should be copy-constructed from that initial State.
  State& operator=(const State& other);
  
  // -----------------------------------------------------------------------------
  // Partial operator=
  //
  // This is hacky, but it works with ATS's multi-state data model to allow
  // subcycling when an entire domain is the thing to be subcycled.
  //
  // This assigns fields and evaluators named domain-*.
  // -----------------------------------------------------------------------------
  void AssignDomain(const State& other, const std::string& domain);

  // Create data structures, finalizing the structure of the state.
  void Setup();

  // Sub-steps in the initialization process. (Used by Amanzi)
  void InitializeEvaluators();
  void InitializeFields();
  void InitializeFieldCopies();
  bool CheckNotEvaluatedFieldsInitialized();
  bool CheckAllFieldsInitialized();

  // Using another state for initialization
  void Initialize(Teuchos::RCP<State> S);
  bool CheckAllFieldsInitialized(Teuchos::RCP<State> S);
  bool CheckNotEvaluatedFieldsInitialized(Teuchos::RCP<State> S);

  // Used by ATS.
  void Initialize();

  // -----------------------------------------------------------------------------
  // State handles mesh management.
  // -----------------------------------------------------------------------------
  // Meshes are "registered" with state.  Creation of meshes is NOT handled by
  // state.
  //
  // Register a mesh under the default key, "domain".
  void RegisterDomainMesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                          bool defoormable=false);

  // Register a mesh under a generic key.
  void RegisterMesh(const Key& key, const Teuchos::RCP<AmanziMesh::Mesh>& mesh,
                    bool deformable=false);

  // Alias a mesh to an existing mesh
  void AliasMesh(const Key& target, const Key& alias);
  bool IsAliasedMesh(const Key& key) const;

  // Remove a mesh.
  void RemoveMesh(const Key& key);

  // Ensure a mesh exists.
  bool HasMesh(const Key& key) const { return GetMesh_(key) != Teuchos::null; }
  bool IsDeformableMesh(const Key& key) const;

  // Mesh accessor.
  Teuchos::RCP<const AmanziMesh::Mesh> GetMesh(const Key& key=Key("domain")) const;
  Teuchos::RCP<AmanziMesh::Mesh> GetDeformableMesh(Key key=Key("domain"));

  // Iterate over meshes.
  typedef MeshMap::const_iterator mesh_iterator;
  mesh_iterator mesh_begin() const { return meshes_.begin(); }
  mesh_iterator mesh_end() const { return meshes_.end(); }
  MeshMap::size_type mesh_count() { return meshes_.size(); }

  // DomainSets are collections of meshes, indexed via NAME_GID and referenced
  // to a parent mesh and sets.
  //
  void RegisterDomainSet(const Key& name,
                         const Teuchos::RCP<AmanziMesh::DomainSet> set);
  bool HasDomainSet(const Key& name) const;
  Teuchos::RCP<const AmanziMesh::DomainSet> GetDomainSet(const Key& name) const;
  

  // -----------------------------------------------------------------------------
  // State handles data management.
  // -----------------------------------------------------------------------------
  // Data is stored and referenced in a common base class, the Field.
  //
  // State manages the creation and consistency of Fields.  Data is "required"
  // of the state.  The require methods act as factories and consistency
  // checks for ownership and type specifiers of the fields.
  //
  // State also manages access to fields.  A Field is "owned" by at most one
  // object -- that object, which is typically either a PK or a
  // FieldEvaluator, may write the solution, and therefore receives non-const
  // pointers to data.  A Field may be used by anyone, but non-owning objects
  // receive const-only pointers to data.  Additionally, fields may be owned
  // by state, meaning that they are independent variables used but not
  // altered by PKs (this is likely changing with the introduction of
  // FieldEvaluators which perform that role).
  //
  // Require Fields from State.
  // -- Require a scalar field, either owned or not.
  void RequireScalar(Key fieldname, Key owner=Key("state"));

  // -- Require a constant vector of given dimension, either owned or not.
  void RequireConstantVector(Key fieldname, Key owner=Key("state"),
                             int dimension=-1);
  void RequireConstantVector(Key fieldname, int dimension=-1);

  // -- Require a vector field, either owned or not.
  Teuchos::RCP<CompositeVectorSpace>
  RequireField(Key fieldname, Key owner="state");

  Teuchos::RCP<CompositeVectorSpace>
  RequireField(Key fieldname, Key owner,
               const std::vector<std::vector<std::string> >& subfield_names);

  // -- A few common, special cases, where we know some of the implied meta-data.
  void RequireGravity();

  // Ensure a mesh exists.
  bool HasField(Key key) const { return GetField_(key) != Teuchos::null; }

  // Field accessor.
  Teuchos::RCP<Field> GetField(Key fieldname, Key pk_name);
  Teuchos::RCP<const Field> GetField(Key fieldname) const;
  void SetField(Key fieldname, Key pk_name, const Teuchos::RCP<Field>& field);

  // Iterate over Fields.
  typedef FieldMap::const_iterator field_iterator;
  field_iterator field_begin() const { return fields_.begin(); }
  field_iterator field_end() const { return fields_.end(); }
  FieldMap::size_type field_count() { return fields_.size(); }

  // Access to Field data
  Teuchos::RCP<const double> GetScalarData(Key fieldname) const;
  Teuchos::RCP<double> GetScalarData(Key fieldname, Key pk_name);

  Teuchos::RCP<const Epetra_Vector> GetConstantVectorData(Key fieldname) const;
  Teuchos::RCP<Epetra_Vector> GetConstantVectorData(Key fieldname, Key pk_name);

  Teuchos::RCP<const CompositeVector> GetFieldData(Key fieldname) const;
  Teuchos::RCP<CompositeVector> GetFieldData(Key fieldname, Key pk_name);

  // Mutator for Field data.
  // -- Modify by pointer, no copy.
  void SetData(Key fieldname, Key pk_name,
                const Teuchos::RCP<double>& data);
  void SetData(Key fieldname, Key pk_name,
                const Teuchos::RCP<Epetra_Vector>& data);
  void SetData(Key fieldname, Key pk_name,
                const Teuchos::RCP<CompositeVector>& data);


  // -----------------------------------------------------------------------------
  // State handles data evaluation.
  // -----------------------------------------------------------------------------
  // To manage lazy yet sufficient updating of models and derivatives of
  // models, we use a graph-based view of data and data dependencies, much
  // like the Phalanx approach.  A directed acyclic graph of dependencies are
  // managed in State, where each node is a FieldEvaluator.
  //
  // -- allows PKs to add to this list to custom evaluators
  Teuchos::ParameterList& FEList() { return state_plist_.sublist("field evaluators"); }
  Teuchos::ParameterList& GetEvaluatorList(const Key& key);
  
  // -- allows PKs to add to this list to initial conditions
  Teuchos::ParameterList& ICList() { return state_plist_.sublist("initial conditions"); }

  // Require FieldEvaluators.
  Teuchos::RCP<FieldEvaluator> RequireFieldEvaluator(Key);
  Teuchos::RCP<FieldEvaluator> RequireFieldEvaluator(Key, Teuchos::ParameterList&);

  // Ensure a FieldEvaluator exists.
  bool HasFieldEvaluator(const Key& key) { return GetFieldEvaluator_(key) != Teuchos::null; }

  // FieldEvaluator accessor.
  Teuchos::RCP<FieldEvaluator> GetFieldEvaluator(const Key&);

  // FieldEvaluator mutator.
  void SetFieldEvaluator(Key key, const Teuchos::RCP<FieldEvaluator>& evaluator);

  // Iterate over evaluators.
  typedef FieldEvaluatorMap::const_iterator evaluator_iterator;
  evaluator_iterator field_evaluator_begin() const { return field_evaluators_.begin(); }
  evaluator_iterator field_evaluator_end() const { return field_evaluators_.end(); }
  FieldEvaluatorMap::size_type field_evaluator_count() { return field_evaluators_.size(); }

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
  // Some fields may have additional copies which corresponds to different times and
  // marked with a timetag
  void RequireTimeTag(Key timetag);
  bool HasTimeTag(Key timetag);
   // Field accessor
  bool HasFieldCopy(Key key, Key timetag);
  Teuchos::RCP<Field> GetFieldCopy(Key fieldname, Key timetag, Key pk_name);
  Teuchos::RCP<const Field> GetFieldCopy(Key fieldname, Key timetag) const;
  void SetFieldCopy(Key fieldname, Key timetag, Key pk_name, const Teuchos::RCP<Field>& field);
  void RequireFieldCopy(Key fieldname,  Key tag, Key copy_owner);
  void CopyField(Key fieldname, Key pk_name, Key timetag);
  Teuchos::RCP<const CompositeVector> GetFieldCopyData(Key fieldname, Key tag) const;
  Teuchos::RCP<CompositeVector> GetFieldCopyData(Key fieldname, Key tag, Key pk_name);

  // Time accessor and mutators.
  double time() const { return time_; }
  void set_time(double new_time);  // note this also evaluates state-owned functions
  void advance_time(double dT) { last_time_ = time(); set_time(time() + dT); }

  double final_time() const { return final_time_; }
  void set_final_time(double new_time) { final_time_ = new_time; }
  double intermediate_time() const { return intermediate_time_; }
  void set_intermediate_time(double new_time) { intermediate_time_ = new_time; }

  double last_time() const { return last_time_; }
  void set_last_time( double last_time) { last_time_ = last_time; }
  double initial_time() const { return initial_time_; }
  void set_initial_time( double initial_time) { initial_time_ = initial_time; }

  // Cycle accessor and mutators.
  int cycle() const { return cycle_; }
  void set_cycle(int cycle) { cycle_ = cycle; }
  void advance_cycle(int dcycle=1) { cycle_ += dcycle; }

  // Position accessor and mutators.
  int position() const { return position_in_tp_; }
  void set_position(int pos ) { position_in_tp_ = pos; }

  // Utility for setting vis flags using blacklist and whitelist
  void InitializeIOFlags();

 private:

  // Accessors that return null if the Key does not exist.
  Teuchos::RCP<AmanziMesh::Mesh> GetMesh_(const Key& key) const;
  Teuchos::RCP<const Field> GetField_(Key fieldname) const;
  Teuchos::RCP<Field> GetField_(Key fieldname);
  Teuchos::RCP<FieldEvaluator> GetFieldEvaluator_(const Key& key);
  Teuchos::RCP<const Functions::MeshPartition> GetMeshPartition_(Key);

  // Consistency checking of fieldnames and types.
  Teuchos::RCP<Field> CheckConsistent_or_die_(Key fieldname,
          FieldType type, Key owner);

  // Containers
  MeshMap meshes_;
  std::map<Key,Key> mesh_aliases_;
  FieldMap fields_;
  FieldFactoryMap field_factories_;
  FieldEvaluatorMap field_evaluators_;

  MeshPartitionMap mesh_partitions_;
  DomainSetMap domain_sets_;

  // meta-data
  double time_;
  double final_time_;
  double intermediate_time_;
  double last_time_;
  double initial_time_;
  std::vector<Key> copy_tag_;


  int cycle_;
  int position_in_tp_;

  // parameter list
  Teuchos::ParameterList state_plist_;
};


// -----------------------------------------------------------------------------
// Non-member functions for I/O of a State.
// -----------------------------------------------------------------------------
// Visualization of State.
void WriteVis(Visualization& vis,
              State& S);

double ReadCheckpoint(State& S,
                      const std::string& filename);

double ReadCheckpointInitialTime(const Comm_ptr_type& comm,
                                 std::string filename);

int ReadCheckpointPosition(const Comm_ptr_type& comm,
                           std::string filename);

void ReadCheckpointObservations(const Comm_ptr_type& comm,
                                std::string filename,
                                Amanzi::ObservationData& obs_data);

void DeformCheckpointMesh(State& S, Key domain);

void WriteStateStatistics(const State& S, const VerboseObject& vo);


}  // namespace Amanzi

#endif
