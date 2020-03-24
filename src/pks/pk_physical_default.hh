/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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

* `"domain name`" ``[string]`` e.g. `"surface`".

  Domains and meshes are 1-to-1, and the empty string refers to the main domain or mesh.  PKs defined on other domains must specify which domain/mesh they refer to.

* `"primary variable key`" ``[string]`` Sets the primary variable.

  The primary variable associated with this PK, i.e. `"pressure`", `"temperature`", `"surface_pressure`", etc.

* `"initial condition`" ``[initial-conditions-spec]``  See InitialConditions_.


   Indicates that the primary variable field has both CELL and FACE objects, and the FACE values are calculated as the average of the neighboring cells.


NOTE: ``PKPhysicalBase (v)-->`` PKDefaultBase_

*/

#ifndef ATS_PK_PHYSICAL_BASE_HH_
#define ATS_PK_PHYSICAL_BASE_HH_



#include "Teuchos_ParameterList.hpp"
#include "TreeVector.hh"

#include "Debugger.hh"

#include "primary_variable_field_evaluator.hh"
#include "PK.hh"
#include "PK_Physical.hh"


namespace Amanzi {

class PK_Physical_Default : public PK_Physical {

  public:
    PK_Physical_Default(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& glist,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~PK_Physical_Default() = default;

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln);
  virtual void Solution_to_State(TreeVector& soln,
                                 const Teuchos::RCP<State>& S);
  virtual void Solution_to_State(const TreeVector& soln,
                                 const Teuchos::RCP<State>& S);


  // new virtual set_states() to also get the primary field evaulator.
  virtual void set_states(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<State>& S_inter,
          const Teuchos::RCP<State>& S_next);

  virtual bool ValidStep();

  // Tag the primary variable as changed in the DAG
  virtual void ChangedSolutionPK(const Teuchos::Ptr<State>& S);
  
  // -- setup
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- initialize
  virtual void Initialize(const Teuchos::Ptr<State>& S);

 protected: // methods

  void DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv);

 protected: // data

  // step validity
  double max_valid_change_;

  // ENORM struct
  typedef struct ENorm_t {
    double value;
    int gid;
  } ENorm_t;
};

} // namespace

#endif
