/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#ifndef ATS_PK_BDF_BASE_HH_
#define ATS_PK_BDF_BASE_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "BDFFnBase.hh"
#include "BDF1_TI.hh"
#include "PK_BDF.hh"
#include "pk_default_base.hh"


namespace Amanzi {

  class PK_BDF_Default : public PK_BDF {

 public:

  // PK_BDF_Default(const Teuchos::RCP<Teuchos::ParameterList>& plist,
  //                Teuchos::ParameterList& FElist,
  //                const Teuchos::RCP<TreeVector>& solution);
  //  :      PKDefaultBase(plist, FElist, solution) {}
  // PK_BDF_Default(Teuchos::ParameterList& FElist,
  //                const Teuchos::RCP<Teuchos::ParameterList>& plist,
  //                const Teuchos::RCP<State>& S,
  //                const Teuchos::RCP<TreeVector>& solution) :
  
  // Virtual destructor
  virtual ~PK_BDF_Default() {}

  // Default implementations of PK methods.
  // -- Setup
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt();

  virtual void set_dt(double dt_);

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  virtual void Solution_to_State(const TreeVector& soln,
                                 const Teuchos::RCP<State>& S);
  virtual void Solution_to_State(TreeVector& soln,
                                 const Teuchos::RCP<State>& S) = 0;

  virtual void set_states(const Teuchos::RCP<const State>& S,
                        const Teuchos::RCP<State>& S_inter,
                        const Teuchos::RCP<State>& S_next);

  // -- ensure a solution is valid
  virtual bool valid_step() { return true; }


  // update the continuation parameter
  virtual void UpdateContinuationParameter(double lambda);

  // -- Check the admissibility of a solution.
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) { return true; }

  // -- Possibly modify the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> up,
          Teuchos::RCP<TreeVector> u) { return false; }

  // -- Possibly modify the correction before it is applied
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

 protected: // data
  // preconditioner assembly control
  bool assemble_preconditioner_;

  // timestep control
  double dt_;
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > time_stepper_;

  // timing
  Teuchos::RCP<Teuchos::Time> step_walltime_;

};

} // namespace

#endif
