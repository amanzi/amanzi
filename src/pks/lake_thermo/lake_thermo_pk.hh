//! Implementation of the PK for the heat equation in the 1D water column (lake)

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#ifndef LAKE_THERMO_PK_HH_
#define LAKE_THERMO_PK_HH_

#include "errors.hh"
#include "pk_bdf_default.hh"
#include "pk_physical_default.hh"

#include "BCs.hh"
#include "Operator.hh"

namespace Amanzi {

class Lake_Thermo_PK : public PK_BDF_Default,
                       public PK_Physical_Default {

 public:

  Lake_Thermo_PK(Teuchos::ParameterList& pk_tree,
                          const Teuchos::RCP<Teuchos::ParameterList>& glist,
                          const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<TreeVector>& solution):
    PK_BDF_Default(pk_tree, glist, S, solution),
    PK_Physical_Default(pk_tree, glist, S, solution),
    PK(pk_tree, glist, S, solution)
  {}


  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next) override;
    
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  virtual void set_dt(double dt) override { dt_ = dt; }

  // initialize.  Note both BDFBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  // Default preconditioner is Picard
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override {
    *Pu = *u;
    return 0;
  }

  // updates the preconditioner, default does nothing
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override {}

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) override;

  virtual bool ValidStep() override {
    return Lake_Thermo_PK::ValidStep() && PK_BDF_Default::ValidStep();
  }
  
  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.

  virtual void ChangedSolution(const Teuchos::Ptr<State>& S) override;

  virtual void ChangedSolution() override;

  virtual double BoundaryValue(const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id);
  virtual int BoundaryDirection(int face_id);
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u);
  
  // PC operator access
  Teuchos::RCP<Operators::Operator> preconditioner() { return preconditioner_; }

  // BC access
  std::vector<int>& bc_markers() { return bc_->bc_model(); }
  std::vector<double>& bc_values() { return bc_->bc_value(); }
  Teuchos::RCP<Operators::BCs> BCs() { return bc_; }

 protected:
  // PC
  Teuchos::RCP<Operators::Operator> preconditioner_;

  // BCs
  Teuchos::RCP<Operators::BCs> bc_;

  // error criteria
  Key conserved_key_;
  Key cell_vol_key_;
  double atol_, rtol_, fluxtol_;

};


} // namespace

#endif
