/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Standard base for most diffusion-dominated PKs, this combines both
domains/meshes of PKPhysicalBase and BDF methods of PKBDFBase.
------------------------------------------------------------------------- */

#ifndef ATS_PK_PHYSICAL_BDF_BASE_HH_
#define ATS_PK_PHYSICAL_BDF_BASE_HH_

#include "errors.hh"
#include "pk_default_base.hh"
#include "pk_bdf_default.hh"
#include "pk_physical_default.hh"

#include "Operator.hh"

namespace Amanzi {

class PK_PhysicalBDF_Default : virtual public PK_BDF_Default, public PK_Physical_Default {

 public:
  PK_PhysicalBDF_Default(Teuchos::ParameterList& FElist,
                          const Teuchos::RCP<Teuchos::ParameterList>& plist,
                          const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<TreeVector>& solution):
    PK_BDF_Default(FElist, plist, S, solution),
    PK_Physical_Default(FElist, plist, S, solution){}


  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);
    
  // virtual void Solution_to_State(TreeVector& solution,
  //                                const Teuchos::RCP<State>& S);
  // virtual void Solution_to_State(const TreeVector& soln,
  //                                const Teuchos::RCP<State>& S);

  virtual void Setup(const Teuchos::Ptr<State>& S);

  virtual std::string name() { return name_; }

  virtual void set_dt(double dt) { dt_ = dt; }

  // initialize.  Note both BDFBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // Default preconditioner is Picard
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
    *Pu = *u;
    return 0;
  }

  // updates the preconditioner, default does nothing
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {}

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution();

  virtual double BoundaryValue(const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id);
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u);

  
  
  // PC operator access
  Teuchos::RCP<Operators::Operator> preconditioner() { return preconditioner_; }

  // BC access
  std::vector<int>& bc_markers() { return bc_markers_; }
  std::vector<double>& bc_values() { return bc_values_; }
  Teuchos::RCP<Operators::BCs> BCs() { return bc_; }

 protected:
  // PC
  Teuchos::RCP<Operators::Operator> preconditioner_;

  // BCs
  std::vector<int> bc_markers_;
  std::vector<double> bc_values_;
  Teuchos::RCP<Operators::BCs> bc_;

  // error criteria
  Key conserved_key_;
  Key cell_vol_key_;
  double atol_, rtol_, fluxtol_;

};


} // namespace

#endif
