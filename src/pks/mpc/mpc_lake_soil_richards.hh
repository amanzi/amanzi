/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/
//! A coupler which integrates temperature models in coupled soil+richards and lake

#ifndef PKS_MPC_LAKE_SOIL_RICHARDS_HH_
#define PKS_MPC_LAKE_SOIL_RICHARDS_HH_

#include "Operator.hh"
#include "pk_physical_bdf_default.hh"
#include "TreeOperator.hh"

#include "strong_mpc.hh"

#include "mpc_coupled_soil.hh"

namespace Amanzi {

class MPCLakeSoilRichards : public StrongMPC<PK_BDF_Default> {
 public:


  MPCLakeSoilRichards(Teuchos::ParameterList& FElist,
				  const Teuchos::RCP<Teuchos::ParameterList>& plist,
				  const Teuchos::RCP<State>& S,
				  const Teuchos::RCP<TreeVector>& soln);

  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
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

  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

 protected:
  // void
  // UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
  //         const Teuchos::RCP<TreeVector>& Pu);

 protected:
  std::string domain_lake_, domain_soil_;

  // sub PKs
  Teuchos::RCP<PK_BDF_Default> lake_pk_;
  Teuchos::RCP<MPCCoupledSoil> soil_mpc_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> lake_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> soil_mesh_;

  // coupled preconditioner
  Teuchos::RCP<Operators::Operator> precon_lake_;
  Teuchos::RCP<Operators::Operator> precon_soil_;
  
  // debugger for dumping vectors
  Teuchos::RCP<Debugger> lake_db_;
  Teuchos::RCP<Debugger> soil_db_;

  Teuchos::RCP<Operators::TreeOperator> op_tree_global_;
  Teuchos::RCP<TreeVector> op_tree_rhs_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCLakeSoilRichards> reg_;

};

} // namespace


#endif
