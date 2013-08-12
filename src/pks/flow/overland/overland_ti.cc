/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"
#include "MatrixMFD_TPFA.hh"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_ICE_FLAG 0
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

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;

    *out_ << "Residual calculation: T0 = " << t_old
          << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << std::setprecision(15);

    S_next_->GetFieldEvaluator("pres_elev")->HasFieldChanged(S_next_.ptr(), name_);
    Teuchos::RCP<const CompositeVector> depth= S_next_->GetFieldData("ponded_depth");
    Teuchos::RCP<const CompositeVector> preselev= S_next_->GetFieldData("pres_elev");

    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  p(" << *c0 <<"): " << (*u)("cell",0,*c0) << std::endl;
        *out_ << "  h(" << *c0 <<"): " << (*depth)("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*depth)("face",fnums0[n]);
        *out_ << std::endl;
        *out_ << "  hz(" << *c0 <<"): " << (*preselev)("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*preselev)("face",fnums0[n]);
        *out_ << std::endl;
        *out_ << " --" << std::endl;
      }
    }
  }
#endif

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    Teuchos::RCP<const CompositeVector> cond =
      S_next_->GetFieldData("upwind_overland_conductivity", name_);
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  cond(" << *c0 << "): " << (*cond)("cell",*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*cond)("face",fnums0[n]);
        *out_ << std::endl;
        *out_ << "  res(" << *c0 << ") (diff): " << (*res)("cell",*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*res)("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif


  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        *out_ << "  res(" << *c0 << ") (acc): " << (*res)("cell",*c0) << std::endl;
      }
    }
  }
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        *out_ << "  res(" << *c0 << ") (source): " << (*res)("cell",*c0) << std::endl;
      }
    }
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
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  p(" << *c0 <<"): " << (*u->data())("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*u->data())("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif

  // apply the preconditioner
  preconditioner_->ApplyInverse(*u, Pu.ptr());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (std::vector<AmanziMesh::Entity_ID>::const_iterator c0=dc_.begin(); c0!=dc_.end(); ++c0) {
      if (*c0 < ncells) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *out_ << "  PC*p(" << *c0 <<"): " << (*Pu->data())("cell",0,*c0);
        for (unsigned int n=0; n!=fnums0.size(); ++n) *out_ << ",  " << (*Pu->data())("face",fnums0[n]);
        *out_ << std::endl;
      }
    }
  }
#endif
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandFlow::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
    *out_ << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  UpdateBoundaryConditionsMarkers_(S_next_.ptr());

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  mfd_preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());
  mfd_preconditioner_->CreateMFDrhsVectors();

  // Patch up BCs in the case of zero conductivity
  FixBCsForPrecon_(S_next_.ptr());

  // 2.b Update local matrices diagonal with the accumulation terms.
  // -- update the accumulation derivatives
  Teuchos::RCP<const CompositeVector> cell_volume =
      S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP<const CompositeVector> depth =
      S_next_->GetFieldData(key_);

  std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = mfd_preconditioner_->Fc_cells();
  unsigned int ncells = cell_volume->size("cell");
  for (unsigned int c=0; c!=ncells; ++c) {
    // accumulation term
    Acc_cells[c] += (*cell_volume)("cell",c) / h;
    Fc_cells[c] += (*depth)("cell",c) * (*cell_volume)("cell",c) / h;
  }

  // Assemble and precompute the Schur complement for inversion.
  // Note boundary conditions are in height variables.
  mfd_preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  if (assemble_preconditioner_) {
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true))
      *out_ << "  assembling..." << std::endl;
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
  */
};


void OverlandFlow::set_preconditioner(const Teuchos::RCP<Operators::Matrix> precon) {
  preconditioner_ = precon;
  mfd_preconditioner_ = Teuchos::rcp_dynamic_cast<Operators::MatrixMFD>(precon);
  ASSERT(mfd_preconditioner_ != Teuchos::null);
  mfd_preconditioner_->set_symmetric(symmetric_);
  mfd_preconditioner_->SymbolicAssembleGlobalMatrices();
  mfd_preconditioner_->CreateMFDmassMatrices(Teuchos::null);
  mfd_preconditioner_->InitPreconditioner();
}

}  // namespace Flow
}  // namespace Amanzi
