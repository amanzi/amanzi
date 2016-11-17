/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! A base class with default implementations of methods for a leaf of the PK tree (a conservation equation, or similar).

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

``PKPhysicalBase`` is a base class providing some functionality for PKs which
are defined on a single mesh, and represent a single process model.  Typically
all leaves of the PK tree will inherit from ``PKPhysicalBase``.

* `"domain`" ``[string]`` **""**, e.g. `"surface`".

  Domains and meshes are 1-to-1, and the empty string refers to the main domain or mesh.  PKs defined on other domains must specify which domain/mesh they refer to.

* `"primary variable`" ``[string]``

  The primary variable associated with this PK, i.e. `"pressure`", `"temperature`", `"surface_pressure`", etc.

* `"initial condition`" ``[initial-condition-spec]``  See InitialConditions_.

  Additionally, the following parameters are supported:

 - `"initialize faces from cell`" ``[bool]`` **false**

   Indicates that the primary variable field has both CELL and FACE objects, and the FACE values are calculated as the average of the neighboring cells.


NOTE: ``PKPhysicalBase (v)-->`` PKDefaultBase_

*/

#ifndef AMANZI_PK_PHYSICAL_BASE_HH_
#define AMANZI_PK_PHYSICAL_BASE_HH_

#include "Debugger.hh"
#include "primary_variable_field_evaluator.hh"
#include "pk_default_base.hh"

namespace Amanzi {

class PKPhysicalBase : public virtual PKDefaultBase {

 public:
  PKPhysicalBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~PKPhysicalBase() {}

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln);
  virtual void solution_to_state(TreeVector& soln,
                                 const Teuchos::RCP<State>& S);

  // new virtual set_states() to also get the primary field evaulator.
  virtual void set_states(const Teuchos::RCP<const State>& S,
          const Teuchos::RCP<State>& S_inter,
          const Teuchos::RCP<State>& S_next);

  virtual bool valid_step();

  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- initialize
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // Accessor for debugger, for use by coupling MPCs
  Teuchos::RCP<Debugger> debugger() { return db_; }

 protected: // methods
  void DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv);

 protected: // data
  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;

  // step validity
  double max_valid_change_;

  // solution and evaluator
  std::string key_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;

  // ENORM struct
  typedef struct ENorm_t {
    double value;
    int gid;
  } ENorm_t;
};

} // namespace

#endif
