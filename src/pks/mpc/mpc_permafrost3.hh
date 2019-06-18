// Rewrite of permafrost PK to simplify
//

#ifndef PKS_MPC_PERMAFROST_THREE_HH_
#define PKS_MPC_PERMAFROST_THREE_HH_

#include "MatrixMFD_Coupled_Surf.hh"
#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA.hh"
#include "mpc_delegate_ewc.hh"
#include "mpc_delegate_water.hh"
#include "pk_physical_bdf_base.hh"

#include "strong_mpc.hh"

namespace Amanzi {

class MPCPermafrost3 : public StrongMPC<PKPhysicalBDFBase> {
 public:

  MPCPermafrost3(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<TreeVector>& soln);

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
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
  void
  UpdateConsistentFaceCorrectionWater_(const Teuchos::Ptr<const TreeVector>& r,
				       const Teuchos::Ptr<const TreeVector>& u,
				       const Teuchos::Ptr<TreeVector>& du);


  int
  ModifyCorrection_FrozenSurface_(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);
  
  // void
  // IteratateFlow_(double h, const Teuchos::RCP<TreeVector>& u);

 protected:
  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_BLOCK_DIAGONAL = 1,
    PRECON_PICARD = 2,
    PRECON_EWC = 3
  };

  // sub PKs
  Teuchos::RCP<PKPhysicalBDFBase> domain_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> domain_energy_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_energy_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // coupled preconditioner
  Teuchos::RCP<TreeMatrix> lin_solver_;
  Teuchos::RCP<Operators::MatrixMFD_Coupled_Surf> precon_;

  // preconditioner methods
  PreconditionerType precon_type_;

  // subblocks of the preconditioner
  // Teuchos::RCP<Operators::MatrixMFD_Surf> pc_flow_;
  // Teuchos::RCP<Operators::MatrixMFD_Surf> pc_energy_;
  Teuchos::RCP<Operators::MatrixMFD> pc_flow_;
  Teuchos::RCP<Operators::MatrixMFD> pc_energy_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> pc_surf_flow_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> pc_surf_energy_;

  // operator for advection in PC
  Teuchos::RCP<Operators::MatrixMFD> pcAdv_;
  Teuchos::RCP<const CompositeVector> adv_field_;
  Teuchos::RCP<const CompositeVector> adv_flux_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWC> surf_ewc_;
  Teuchos::RCP<MPCDelegateEWC> sub_ewc_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;
  bool consistent_cells_;
  

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost3> reg_;

};

} // namespace


#endif
