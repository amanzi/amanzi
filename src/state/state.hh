/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for the State.  State is a simple data-manager, allowing PKs to
   require, read, and write various fields, including:
    -- Acts as a factory for fields through the various require methods.
    -- Provides some data protection by providing both const and non-const
       data pointers to PKs.
    -- Provides some initialization capability -- this is where all
       independent variables can be initialized (as independent variables
       are owned by state, not by any PK).
   ------------------------------------------------------------------------- */

#ifndef STATE_STATE_HH_
#define STATE_STATE_HH_

#include <string>
#include <vector>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"

#include "composite_vector.hh"
#include "composite_vector_factory.hh"

#include "state_defs.hh"

#include "visualization.hh"
#include "checkpoint.hh"

#include "field.hh"
#include "field_scalar.hh"
#include "field_constant_vector.hh"
#include "field_composite_vector.hh"

namespace Amanzi {

class FieldEvaluator;

class State {

 private:

  typedef std::map<Key, Teuchos::RCP<const AmanziMesh::Mesh> > MeshMap;
  typedef std::map<Key, Teuchos::RCP<CompositeVectorFactory> > FieldFactoryMap;
  typedef std::map<Key, Teuchos::RCP<Field> > FieldMap;
  typedef std::map<Key, Teuchos::RCP<FieldEvaluator> > FieldEvaluatorMap;

 public:

  // Default constructor.
  State();

  // Usual constructor.
  explicit State(Teuchos::ParameterList& state_plist);

  // Copy constructor, copies memory not pointers.
  explicit State(const State& other);

  // Assignment operator, copies memory not pointers.  Note this
  // implementation requires the State being copied has the same structure (in
  // terms of fields, order of fields, etc) as *this.  This really means that
  // it should be a previously-copy-constructed version of the State.  One and
  // only one State should be instantiated and populated -- all other States
  // should be copy-constructed from that initial State.
  State& operator=(const State& other);

  // Create data structures, finalizing the structure of the state.
  void Setup();

  // Initialize field evaluators using ICs set by PKs.
  void Initialize();


  // -----------------------------------------------------------------------------
  // State handles mesh management.
  // -----------------------------------------------------------------------------
  // Meshes are "registered" with state.  Creation of meshes is NOT handled by
  // state.
  //
  // Register a mesh under the default key, "domain".
  void RegisterDomainMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // Register a mesh under a generic key.
  void RegisterMesh(Key key, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // Remove a mesh.
  void RemoveMesh(Key key);

  // Ensure a mesh exists.
  bool HasMesh(Key key) const { return GetMesh_(key) != Teuchos::null; }

  // Mesh accessor.
  Teuchos::RCP<const AmanziMesh::Mesh> GetMesh(Key key=Key("domain")) const;

  // Iterate over meshes.
  typedef MeshMap::const_iterator mesh_iterator;
  mesh_iterator mesh_begin() const { return meshes_.begin(); }
  mesh_iterator mesh_end() const { return meshes_.end(); }
  MeshMap::size_type mesh_count() { return meshes_.size(); }

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
  Teuchos::RCP<CompositeVectorFactory>
  RequireField(Key fieldname, Key owner="state");

  // -- A few common, special cases, where we know some of the implied meta-data.
  void RequireGravity();

  // Ensure a mesh exists.
  bool HasField(Key key) const { return GetField_(key) != Teuchos::null; }

  // Field accessor.
  Teuchos::RCP<Field> GetField(Key fieldname, Key pk_name);
  Teuchos::RCP<const Field> GetField(Key fieldname) const;

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
  // Require FieldEvaluators.
  Teuchos::RCP<FieldEvaluator> RequireFieldEvaluator(Key);
  Teuchos::RCP<FieldEvaluator> RequireFieldEvaluator(Key, Teuchos::ParameterList&);

  // Ensure a FieldEvaluator exists.
  bool HasFieldEvaluator(Key key) { return GetFieldEvaluator_(key) != Teuchos::null; }

  // FieldEvaluator accessor.
  Teuchos::RCP<FieldEvaluator> GetFieldEvaluator(Key);

  // FieldEvaluator mutator.
  void SetFieldEvaluator(Key key, const Teuchos::RCP<FieldEvaluator>& evaluator);

  // Iterate over evaluators.
  typedef FieldEvaluatorMap::const_iterator evaluator_iterator;
  evaluator_iterator field_evaluator_begin() const { return field_evaluators_.begin(); }
  evaluator_iterator field_evaluator_end() const { return field_evaluators_.end(); }
  FieldEvaluatorMap::size_type field_evaluator_count() { return field_evaluators_.size(); }


  // -----------------------------------------------------------------------------
  // State handles model parameters.
  // -----------------------------------------------------------------------------
  // Some model parameters may be common to many PKs, Evaluators, boundary
  // conditions, etc.  Access to the parameters required to make these models
  // is handled through state.
  //
  // Get a parameter list.
  Teuchos::ParameterList GetModelParameters(std::string modelname); 


  // -----------------------------------------------------------------------------
  // State is representative of an instant in time and a single cycle within
  // the time integration process.
  // -----------------------------------------------------------------------------
  // Time accessor and mutators.
  double time() const { return time_; }
  void set_time(double new_time);  // note this also evaluates state-owned functions
  void advance_time(double dT) { set_time(time() + dT); }

  // Cycle accessor and mutators.
  int cycle() const { return cycle_; }
  void set_cycle(int cycle) { cycle_ = cycle; }
  void advance_cycle(int dcycle=1) { cycle_ += dcycle; }

private:
  // sub-steps in the initialization process.
  void InitializeEvaluators_();
  void InitializeFields_();
  bool CheckNotEvaluatedFieldsInitialized_();
  bool CheckAllFieldsInitialized_();

  // Accessors that return null if the Key does not exist.
  Teuchos::RCP<const AmanziMesh::Mesh> GetMesh_(Key key) const;
  Teuchos::RCP<const Field> GetField_(Key fieldname) const;
  Teuchos::RCP<Field> GetField_(Key fieldname);
  Teuchos::RCP<FieldEvaluator> GetFieldEvaluator_(Key key);

  // Consistency checking of fieldnames and types.
  Teuchos::RCP<Field> CheckConsistent_or_die_(Key fieldname,
          FieldType type, Key owner);

  // Containers
  MeshMap meshes_;
  FieldMap fields_;
  FieldFactoryMap field_factories_;
  FieldEvaluatorMap field_evaluators_;

  // meta-data
  double time_;
  int cycle_;

  // parameter list
  Teuchos::ParameterList state_plist_;
};


// -----------------------------------------------------------------------------
// Non-member functions for I/O of a State.
// -----------------------------------------------------------------------------
// Visualization of State.
void WriteVis(const Teuchos::Ptr<Visualization>& vis,
              const Teuchos::Ptr<State>& S);

// Checkpointing State.
void WriteCheckpoint(const Teuchos::Ptr<Checkpoint>& ckp,
                     const Teuchos::Ptr<State>& S,
                     double dt);

double ReadCheckpoint(Epetra_MpiComm* comm,
                      const Teuchos::Ptr<State>& S,
                      std::string filename);

double ReadCheckpointInitialTime(Epetra_MpiComm* comm,
                      std::string filename);

void DeformCheckpointMesh(const Teuchos::Ptr<State>& S);

} // namespace amanzi

#endif
