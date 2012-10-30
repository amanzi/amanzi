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

#define DEBUG_FLAG 0

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
  S_inter_->set_time(t_old);
  S_next_ ->set_time(t_new);
  double h = t_new - t_old;
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
  ApplyDiffusion_(S_next_, res);

#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after diffusion): " << (*res)("cell",0,0) << " "
              << (*res)("face",0,0) << std::endl;
  }
#endif

  // accumulation term
  AddAccumulation_(res);
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  res0 (after accumulation): " << (*res)("cell",0,0) << " "
              << (*res)("face",0,0) << std::endl;
  }
#endif

  // add rhs load value
  AddLoadValue_(res);
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

  // check for nans in residual
  Teuchos::RCP<const CompositeVector> res = u->data();
  for (int c=0; c!=res->size("cell"); ++c) {
    if (boost::math::isnan<double>((*res)("cell",c))) {
      int mypid = S_next_->GetMesh("surface")->get_comm()->MyPID();
      *out_ << "Cutting time step due to NaN in cell residual on proc "
                << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=res->size("face"); ++f) {
    if (boost::math::isnan<double>((*res)("face",f))) {
      int mypid = S_next_->GetMesh("surface")->get_comm()->MyPID();
      *out_ << "Cutting time step due to NaN in face residual on proc "
                << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  // Dump correction
#if DEBUG_FLAG
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " "
              << (*Pu->data())("face",0,0) << std::endl;
  }
#endif

  // check for nans in preconditioned residual
  Teuchos::RCP<const CompositeVector> Pres = Pu->data();
  for (int c=0; c!=Pres->size("cell"); ++c) {
    if (boost::math::isnan<double>((*Pres)("cell",c))) {
      int mypid = S_next_->GetMesh("surface")->get_comm()->MyPID();
      *out_ << "Cutting time step due to NaN in PC'd cell residual on proc "
                << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=Pres->size("face"); ++f) {
    if (boost::math::isnan<double>((*Pres)("face",f))) {
      int mypid = S_next_->GetMesh("surface")->get_comm()->MyPID();
      *out_ << "Cutting time step due to NaN in PC'd face residual on proc "
                << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
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
  S_next_->set_time(t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_    ->Compute(S_next_->time());
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
  preconditioner_->CreateMFDstiffnessMatrices(cond.ptr());
  preconditioner_->CreateMFDrhsVectors();

  // 2. Update local matrices diagonal with the accumulation terms.
  Teuchos::RCP<const CompositeVector> cell_volume =
      S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP<const CompositeVector> pres =
      S_next_->GetFieldData(key_);

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  int ncells = cell_volume->size("cell");
  for (int c=0; c!=ncells; ++c) {
    // accumulation term
    Acc_cells[c] += (*cell_volume)("cell",c) / h;
    Fc_cells[c] += (*pres)("cell",c) * (*cell_volume)("cell",c) / h;
  }

  // Assemble and precompute the Schur complement for inversion.
  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();
  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);

  // dump the schur complement
  //  Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
  //  std::stringstream filename_s;
  //  filename_s << "schur_" << S_next_->cycle() << ".txt";
  //  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  //  *out_ << "updated precon " << S_next_->cycle() << std::endl;


  preconditioner_->UpdatePreconditioner();

  /*
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

    preconditioner_->Apply(*dp->data(), df->data());
    double df_loc = (*df->data())("cell",c);
    df->Update(-1.0, *f2, 1.0, *f1, 1.0);
    double error = enorm(f1, df);
    //
      AmanziGeometry::Point point = S_next_->GetMesh("surface")->cell_centroid(c);
      if (error > 1e-5) {
        std::cout << "Bad error at cell: " << c << std::endl;
        AmanziMesh::Entity_ID_List faces;
        std::vector<int> fdirs;
        S_next_->GetMesh("surface")->cell_get_faces_and_dirs(c, &faces, &fdirs);
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
