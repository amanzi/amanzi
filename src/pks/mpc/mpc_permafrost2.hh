// Rewrite of permafrost PK to simplify
//

#ifndef PKS_MPC_PERMAFROST2_HH_
#define PKS_MPC_PERMAFROST2_HH_

#include "MatrixMFD_Coupled_Surf.hh"
#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA.hh"
#include "mpc_delegate_ewc.hh"
#include "mpc_delegate_water.hh"
#include "pk_physical_bdf_base.hh"

#include "strong_mpc.hh"

namespace Amanzi {

class MPCPermafrost2 : public StrongMPC<PKPhysicalBDFBase> {
 public:

  MPCPermafrost2(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<TreeVector>& soln);

  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk fun().
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Modify the predictor.
  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual bool modify_correction(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

 protected:
  void
  UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
          const Teuchos::RCP<TreeVector>& Pu);


  int
  ModifyCorrection_FrozenSurface_(double h, Teuchos::RCP<const TreeVector> res,
          Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);
  
  // void
  // IteratateFlow_(double h, const Teuchos::RCP<TreeVector>& u);

 protected:

  // sub PKs
  Teuchos::RCP<PKPhysicalBDFBase> domain_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> domain_energy_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_energy_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // coupled preconditioner
  Teuchos::RCP<Operators::MatrixMFD_Coupled_Surf> precon_;
  Teuchos::RCP<Operators::MatrixMFD_Surf> pc_flow_;
  Teuchos::RCP<Operators::MatrixMFD_Surf> pc_energy_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> pc_surf_flow_;
  Teuchos::RCP<Operators::MatrixMFD_TPFA> pc_surf_energy_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWC> ewc_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;
  bool consistent_cells_;
  

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost2> reg_;

};

} // namespace


#endif
