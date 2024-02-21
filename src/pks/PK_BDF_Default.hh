/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon
*/

/*
  Default implementation helper class for implementing PKs that are also BDF
  time integrable.
*/

#pragma once

#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "Key.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "PDE_HelperDiscretization.hh"
#include "EvaluatorPrimary.hh"

#include "PK_Physical.hh"

namespace Amanzi {

class PK_BDF_Default : public PK_BDF {
 public:

  PK_BDF_Default(const Comm_ptr_type& comm,
                 Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S);

  // access to operators and PDEs in sub-PKs
  virtual Teuchos::RCP<Operators::Operator> getOperator(const Operators::Operator_kind& type) override {
    return Teuchos::null;
  }

  // -- Choose a time step compatible with physics.
  virtual double getDt() override;
  virtual void setDt(double dt) override;

  // -- PK API
  virtual void parseParameterList() override;
  virtual void setup() override;
  virtual void initialize() override;
  virtual bool advanceStep(double t_old, double t_new, bool reinit) override;
  virtual void commitStep(double t_old, double t_new, const Tag& tag) override;

  // -- BDF API
  virtual bool isAdmissible(Teuchos::RCP<const TreeVector> up) override { return true; }
  virtual bool modifyPredictor(double h, Teuchos::RCP<const TreeVector> up, Teuchos::RCP<TreeVector> u) override {
    return false;
  }

  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  modifyCorrection(double h,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override
  {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  virtual void markChangedSolution() override {
    markChangedSolution(tag_next_);
  }

 protected:                      // data
  bool assemble_preconditioner_; // preconditioner assembly control
  bool strongly_coupled_;        // if we are coupled, no need to make a TI

  // timestep control, solution vector
  Teuchos::RCP<TreeVector> solution_;
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> time_stepper_;

  // timing
  Teuchos::RCP<Teuchos::Time> step_walltime_;
};

} // namespace Amanzi

