/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_RES_FLAG 0


// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void OverlandFlow::fun( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g ) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  // bookkeeping
  double h = t_new - t_old;
  ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  Teuchos::RCP<CompositeVector> u = u_new->data();

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "OverlandFlow Residual calculation:" << std::endl;
    *out_ << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,0) << std::endl;
  }
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after diffusion): " << (*res)("cell",0,0) << " "
              << (*res)("face",0,0) << std::endl;
  }
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after accumulation): " << (*res)("cell",0,0) << " "
              << (*res)("face",0,0) << std::endl;
  }
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after source): " << (*res)("cell",0,0) << " "
              << (*res)("face",0,0) << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void OverlandFlow::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,0) << " "
              << (*u->data())("face",0,0) << std::endl;
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u, Pu.ptr());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " "
              << (*Pu->data())("face",0,0) << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandFlow::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon update at t = " << t << std::endl;
    *out_ << "  p0: " << (*up->data())("cell",0,0) << " "
              << (*up->data())("face",0,0) << std::endl;
  }
#endif

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditionsNoElev_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  conductivity0: " << (*cond)("cell",0,0) << " "
              << (*cond)("face",0,0) << std::endl;
  }
#endif

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  mfd_preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());
  mfd_preconditioner_->CreateMFDrhsVectors();

  // 2. Update local matrices diagonal with the accumulation terms.
  Teuchos::RCP<const CompositeVector> cell_volume =
      S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP<const CompositeVector> pres =
      S_next_->GetFieldData(key_);

  std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = mfd_preconditioner_->Fc_cells();
  int ncells = cell_volume->size("cell");
  for (int c=0; c!=ncells; ++c) {
    // accumulation term
    Acc_cells[c] += (*cell_volume)("cell",c) / h;
    Fc_cells[c] += (*pres)("cell",c) * (*cell_volume)("cell",c) / h;
  }

  // Assemble and precompute the Schur complement for inversion.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (assemble_preconditioner_) {
    mfd_preconditioner_->AssembleGlobalMatrices();
    mfd_preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);
    mfd_preconditioner_->UpdatePreconditioner();
  }

  /*
  // dump the schur complement
  Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  *out_ << "updated precon " << S_next_->cycle() << std::endl;


  // print the rel perm
  Teuchos::RCP<const CompositeVector> num_rel_perm =
      S_next_->GetFieldData("upwind_overland_conductivity");
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("overland_conductivity");
  *out_ << "REL PERM: " << std::endl;
  rel_perm->Print(*out_);
  *out_ << std::endl;
  *out_ << "UPWINDED REL PERM: " << std::endl;
  num_rel_perm->Print(*out_);
  */
};


void OverlandFlow::set_preconditioner(const Teuchos::RCP<Operators::Matrix> precon) {
  preconditioner_ = precon;
  mfd_preconditioner_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_preconditioner_ != Teuchos::null);
  mfd_preconditioner_->SetSymmetryProperty(symmetric_);
  mfd_preconditioner_->SymbolicAssembleGlobalMatrices();
  mfd_preconditioner_->CreateMFDmassMatrices(Teuchos::null);
  mfd_preconditioner_->InitPreconditioner();
}


// Runs a very expensive FD test of the Jacobian and prints out an enorm
//  measure of the error.
void OverlandFlow::test_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  Teuchos::RCP<TreeVector> dp = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> df = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> uold = Teuchos::rcp(new TreeVector(*up));
  Teuchos::RCP<TreeVector> unew = Teuchos::rcp(new TreeVector(*up));

  double maxval = 0.0;

  std::cout.precision(15);

  int ncells = up->data()->size("cell");
  for (int c=0; c!=ncells; ++c) {
    *unew = (*up);
    fun(t-h, t, uold, unew, f1);

    dp->PutScalar(0.0);
    (*dp->data())("cell",c) = 0.00001;
    unew->Update(1.0, *dp, 1.0);
    fun(t-h, t, uold, unew, f2);

    preconditioner_->Apply(*dp, df.ptr());
    double df_loc = (*df->data())("cell",c);
    df->Update(-1.0, *f2, 1.0, *f1, 1.0);
    double error = enorm(f1, df);
    //
      AmanziGeometry::Point point = mesh_->cell_centroid(c);
      if (error > 1e-5) {
        std::cout << "Bad error at cell: " << c << std::endl;
        AmanziMesh::Entity_ID_List faces;
        std::vector<int> fdirs;
        mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
        std::cout << "faces: " << faces[0] << ", "
            << faces[1] << ", "
            << faces[2] << ", "
            << faces[3] << std::endl;

      }

      std::cout << "error: " << std::scientific << (*df->data())("cell",c) << std::endl;
      std::cout << "  cell center: " << point << std::endl;
      std::cout << "  f_1: " << std::scientific << (*f1->data())("cell",c) << std::endl;
      std::cout << "  f_2: " << std::scientific << (*f2->data())("cell",c) << std::endl;
      std::cout << "  df:  " << std::scientific << df_loc << std::endl;
      //    }
    maxval = std::max(maxval, error);
  }
  std::cout << "Testing PC with FD.  Error: " << maxval << std::endl;
};

}  // namespace Flow
}  // namespace Amanzi
