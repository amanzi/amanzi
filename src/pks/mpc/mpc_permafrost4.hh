// Rewrite of permafrost PK, new operators, based on mpc_subsurface
//

#ifndef PKS_MPC_PERMAFROST_FOUR_HH_
#define PKS_MPC_PERMAFROST_FOUR_HH_

#include "mpc_delegate_ewc.hh"
#include "mpc_delegate_water.hh"
#include "mpc_subsurface.hh"

namespace Amanzi {

class MPCPermafrost4 : public MPCSubsurface {
 public:

  MPCPermafrost4(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<TreeVector>& soln);

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk Functional().
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- Apply preconditioner to r and returns the result in Pr.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r, Teuchos::RCP<TreeVector> Pr);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> r,
                       Teuchos::RCP<const TreeVector> u, 
                       Teuchos::RCP<TreeVector> du);

 protected:
  // sub PKs
  Teuchos::RCP<PKPhysicalBDFBase> domain_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> domain_energy_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_energy_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // off-diagonal terms
  Teuchos::RCP<Operators::OperatorAccumulation> dE_dp_surf_;

  
  // EWC delegate for the surface
  //  Teuchos::RCP<MPCDelegateEWC> surf_ewc_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost4> reg_;
  Key domain_surf, domain_ss;

};

} // namespace


#endif
