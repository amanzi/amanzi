// Rewrite of permafrost PK to simplify
//

#ifndef PKS_MPC_COUPLED_WATER_HH_
#define PKS_MPC_COUPLED_WATER_HH_

#include "Operator.hh"
#include "mpc_delegate_water.hh"
#include "pk_physical_bdf_base.hh"

#include "strong_mpc.hh"

namespace Amanzi {

class MPCCoupledWater : public StrongMPC<PKPhysicalBDFBase> {
 public:

  MPCCoupledWater(const Teuchos::RCP<Teuchos::ParameterList>& plist,
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

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

 protected:
  // void
  // UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
  //         const Teuchos::RCP<TreeVector>& Pu);

 protected:

  // sub PKs
  Teuchos::RCP<PKPhysicalBDFBase> domain_flow_pk_;
  Teuchos::RCP<PKPhysicalBDFBase> surf_flow_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // coupled preconditioner
  Teuchos::RCP<Operators::Operator> precon_;
  Teuchos::RCP<Operators::Operator> precon_surf_;
  Teuchos::RCP<Operators::Operator> lin_solver_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;
  bool consistent_cells_;
  
  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory_ATS<MPCCoupledWater> reg_;

};

} // namespace


#endif
