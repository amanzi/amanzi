/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
BDF.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_BDF_BASE_HH_
#define AMANZI_PK_BDF_BASE_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "BDFFnBase.hh"
#include "BDF1_TI.hh"
#include "Matrix.hh"
#include "pk_default_base.hh"

namespace Amanzi {

class PKBDFBase : public virtual PKDefaultBase,
                  public BDFFnBase<TreeVector> {

 public:

  PKBDFBase(Teuchos::ParameterList& plist,
            const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution),
      backtracking_(false) {}

  // Virtual destructor
  virtual ~PKBDFBase() {}

  // Default implementations of PK methods.
  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- initialize
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt();

  // Operator access/mutate/use
  virtual Teuchos::RCP<Operators::Matrix> preconditioner() { return preconditioner_; }
  virtual void set_preconditioner(const Teuchos::RCP<Operators::Matrix> preconditioner) {
    preconditioner_ = preconditioner; }
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt);

  // -- Check the admissibility of a solution.
  virtual bool is_admissible(Teuchos::RCP<const TreeVector> up) { return true; }

  // -- Possibly modify the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator.
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up) { return false; }

  // -- Possibly modify the correction before it is applied
  virtual bool modify_correction(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
    return false; }

 protected: // data
  // timestep control
  double dt_;
  Teuchos::RCP<BDF1_TI<TreeVector> > time_stepper_;
  bool backtracking_;
  int backtracking_count_;
  int backtracking_iterations_;

  // operators
  Teuchos::RCP<Operators::Matrix> preconditioner_;
  bool assemble_preconditioner_;

  // timing
  Teuchos::RCP<Teuchos::Time> step_walltime_;

};

} // namespace

#endif
