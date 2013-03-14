/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for a StrongMPC which uses a preconditioner in which the
   block-diagonal cell-local matrix is dense.  If the system looks something
   like:

   A( y1, y2, x, t ) = 0
   B( y1, y2, x, t ) = 0

   where y1,y2 are spatially varying unknowns that are discretized using the MFD
   method (and therefore have both cell and face unknowns), an approximation to
   the Jacobian is written as

   [  dA_c/dy1_c  dA_c/dy1_f   dA_c/dy2_c       0      ]
   [  dA_f/dy1_c  dA_f/dy1_f      0              0      ]
   [  dB_c/dy1_c     0          dB_c/dy2_c  dB_c/dy2_f ]
   [      0           0          dB_f/dy2_c  dB_f/dy2_f ]


   Note that the upper left block is the standard preconditioner for the A
   system, and the lower right block is the standard precon for the B system,
   and we have simply added cell-based couplings, dA_c/dy2_c and dB_c/dy1_c.

   In the temperature/pressure system, these correspond to d_water_content /
   d_temperature and d_energy / d_pressure.

   ------------------------------------------------------------------------- */

#include <fstream>

#include "field_evaluator.hh"
#include "matrix_mfd.hh"
#include "matrix_coupled_mfd.hh"
#include "mpc_coupled_cells.hh"

namespace Amanzi {

MPCCoupledCells::MPCCoupledCells(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist,soln),
    StrongMPC(plist,soln) {}

void MPCCoupledCells::setup(const Teuchos::Ptr<State>& S) {
  StrongMPC::setup(S);

  A_key_ = plist_.get<std::string>("conserved quantity A");
  B_key_ = plist_.get<std::string>("conserved quantity B");
  y1_key_ = plist_.get<std::string>("primary variable A");
  y2_key_ = plist_.get<std::string>("primary variable B");
  dA_dy2_key_ = std::string("d")+A_key_+std::string("_d")+y2_key_;
  dB_dy1_key_ = std::string("d")+B_key_+std::string("_d")+y1_key_;

  Key mesh_key = plist_.get<std::string>("mesh key");
  mesh_ = S->GetMesh(mesh_key);

  // Create the precon
  Teuchos::ParameterList pc_sublist = plist_.sublist("Coupled PC");
  mfd_preconditioner_ =
      Teuchos::rcp(new Operators::MatrixCoupledMFD(pc_sublist, mesh_));

  // Set the sub-blocks from the sub-PK's preconditioners.
  Teuchos::RCP<Operators::Matrix> pcA = sub_pks_[0]->preconditioner();
  Teuchos::RCP<Operators::Matrix> pcB = sub_pks_[1]->preconditioner();

#ifdef ENABLE_DBC
  Teuchos::RCP<Operators::MatrixMFD> pcA_mfd =
      Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(pcA);
  ASSERT(pcA_mfd != Teuchos::null);
  Teuchos::RCP<Operators::MatrixMFD> pcB_mfd =
      Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(pcB);
  ASSERT(pcB_mfd != Teuchos::null);
#else
  Teuchos::RCP<Operators::MatrixMFD> pcA_mfd =
      Teuchos::rcp_static_cast<Operators::MatrixMFD>(pcA);
  Teuchos::RCP<Operators::MatrixMFD> pcB_mfd =
      Teuchos::rcp_static_cast<Operators::MatrixMFD>(pcB);
#endif

  mfd_preconditioner_->SetSubBlocks(pcA_mfd, pcB_mfd);

  // setup and initialize the preconditioner
  mfd_preconditioner_->SymbolicAssembleGlobalMatrices();
  mfd_preconditioner_->InitPreconditioner();
  preconditioner_ = mfd_preconditioner_;
}


// updates the preconditioner
void MPCCoupledCells::update_precon(double t, Teuchos::RCP<const TreeVector> up,
        double h) {
  StrongMPC::update_precon(t,up,h);

  // Update and get the off-diagonal terms.
  S_next_->GetFieldEvaluator(A_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y2_key_);
  S_next_->GetFieldEvaluator(B_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, y1_key_);
  Teuchos::RCP<const CompositeVector> dA_dy2 = S_next_->GetFieldData(dA_dy2_key_);
  Teuchos::RCP<const CompositeVector> dB_dy1 = S_next_->GetFieldData(dB_dy1_key_);

  // scale by 1/h
  Epetra_MultiVector Ccc(*dA_dy2->ViewComponent("cell",false));
  Ccc = *dA_dy2->ViewComponent("cell",false);
  Ccc.Scale(1./h);

  Epetra_MultiVector Dcc(*dB_dy1->ViewComponent("cell",false));
  Dcc = *dB_dy1->ViewComponent("cell",false);
  Dcc.Scale(1./h);

  // Assemble the precon, form Schur complement
  mfd_preconditioner_->ComputeSchurComplement(Ccc, Dcc);
  mfd_preconditioner_->UpdatePreconditioner();

}


// applies preconditioner to u and returns the result in Pu
void MPCCoupledCells::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // std::ofstream file;
  // file.open("residual.txt");
  // u->Print(file);

  preconditioner_->ApplyInverse(*u, Pu.ptr());

  // Pu->Print(file);

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    AmanziMesh::Entity_ID_List fnums,fnums0;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(c0_, &fnums0, &dirs);
    mesh_->cell_get_faces_and_dirs(c1_, &fnums, &dirs);

    Teuchos::OSTab tab = getOSTab();
    *out_ << "Preconditioned Updates:" << std::endl;
    *out_ << "  Pp0: " << (*Pu->SubVector(0)->data())("cell",c0_) << " " << (*Pu->SubVector(0)->data())("face",fnums0[0]) << std::endl;
    *out_ << "  Pp1: " << (*Pu->SubVector(0)->data())("cell",c1_) << " " << (*Pu->SubVector(0)->data())("face",fnums[1]) << std::endl;
    *out_ << "  PT0: " << (*Pu->SubVector(1)->data())("cell",c0_) << " " << (*Pu->SubVector(1)->data())("face",fnums0[0]) << std::endl;
    *out_ << "  PT1: " << (*Pu->SubVector(1)->data())("cell",c1_) << " " << (*Pu->SubVector(1)->data())("face",fnums[1]) << std::endl;
  }
}


} //  namespace
