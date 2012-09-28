/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for flow and energy.  This couples using a
block-diagonal coupler.
------------------------------------------------------------------------- */

#ifndef MPC_PREC_COUPLED_FLOW_ENERGY_HH_
#define MPC_PREC_COUPLED_FLOW_ENERGY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "state.hh"
#include "strong_mpc.hh"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVbrMatrix.h"

#include "tree_vector.hh"
#include "bdf_fn_base.hh"
#include "bdf_time_integrator.hh"

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

#include "richards.hh"
#include "two_phase.hh"

#include "ml_MultiLevelPreconditioner.h"

#include "Ifpack.h"
#include "Ifpack_ILU.h"
#include "Ifpack_AdditiveSchwarz.h"


// class Amanzi::Flow::Richards;

namespace Amanzi {

class MPCCoupledFlowEnergy : public StrongMPC {

 public:
  MPCCoupledFlowEnergy(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, soln),
      StrongMPC(plist, soln) {}

  // initialize the preconditioner
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);
  
  void ComputeShurComplementPK();

 protected:
//   void precon_version1(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  void SymbolicAssembleGlobalMatrices(const Teuchos::Ptr<State>& S);
  
  Teuchos::RCP<Epetra_MultiVector> D_pT_;
  Teuchos::RCP<Epetra_MultiVector> D_Tp_;

//   PreconMethod method_;
  double damping_;

 private:
         
  void InitPreconditioner(Teuchos::ParameterList& prec_plist);
         
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  Teuchos::RCP<Epetra_FEVbrMatrix> P2f2f_;
  Teuchos::RCP<Epetra_VbrMatrix> A2f2p_;
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Cell_Couple_Inv_ ;
  Teuchos::RCP<const Epetra_Map> supermap_;
  
  Teuchos::RCP<Amanzi::Flow::Richards> flow_pk;
  Teuchos::RCP<Amanzi::Energy::TwoPhase> energy_pk;
  
  Teuchos::RCP<const Epetra_BlockMap> fmap_wghost;
  Teuchos::RCP<const Epetra_BlockMap> double_fmap;
  Teuchos::RCP<const Epetra_BlockMap> double_cmap;
  Teuchos::RCP<const Epetra_BlockMap> double_fmap_wghost;
  
  bool is_matrix_constructed;
  bool decoupled;
  
   enum { TRILINOS_ML, TRILINOS_ILU, TRILINOS_BLOCK_ILU, HYPRE_AMG, HYPRE_EUCLID, HYPRE_PARASAILS } prec_method_; 
  
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec_;
  Teuchos::ParameterList ml_plist_, coupled_pc_;

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff_;
  Teuchos::ParameterList hypre_plist_;
  int hypre_ncycles_, hypre_nsmooth_;
  double hypre_tol_, hypre_strong_threshold_;
#endif

  Teuchos::RCP<Ifpack_ILU> ilu_prec_;
  Teuchos::ParameterList ilu_plist_;

  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;
  
  // factory registration
  static RegisteredPKFactory<MPCCoupledFlowEnergy> reg_;

};

} // namespace


#endif
