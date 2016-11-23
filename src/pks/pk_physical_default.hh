/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a physical PK.
------------------------------------------------------------------------- */

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
  virtual ~PK_Physical_Default() {}

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln);
  virtual void Solution_to_State(TreeVector& soln,
                                 const Teuchos::RCP<State>& S);
  virtual void Solution_to_State(const TreeVector& soln,
                                 const Teuchos::RCP<State>& S);


  // new virtual set_states() to also get the primary field evaulator.
  virtual void set_states(const Teuchos::RCP<const State>& S,
          const Teuchos::RCP<State>& S_inter,
          const Teuchos::RCP<State>& S_next);

  virtual bool ValidStep();

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
