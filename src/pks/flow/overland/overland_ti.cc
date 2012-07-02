/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Gianmarco Manzini
         Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "overland.hh"

namespace Amanzi {
namespace Flow {

#if 0
}}
#endif

#define debug_flag 0


// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void OverlandFlow::fun( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g ) {
  S_inter_->set_time(t_old);
  S_next_ ->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
#if debug_flag
  std::cout << "OverlandFlow Residual calculation:" << std::endl;
  std::cout << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u)("cell",0,9) << " " << (*u)("face",0,29) << std::endl;
#endif
  //print_vector(S_next_,u,"fun") ;

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_);

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_, res);
#if debug_flag
  std::cout << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after diffusion): " << (*res)("cell",0,9) << " " << (*res)("face",0,29) << std::endl;
#endif

  // accumulation term
  AddAccumulation_(res);
#if debug_flag
  std::cout << "  res0 (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after accumulation): " << (*res)("cell",0,9) << " " << (*res)("face",0,29) << std::endl;
#endif

  // add rhs load value
  AddLoadValue_(S_next_,res);
#if debug_flag
  std::cout << "  res0 (after source): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
  std::cout << "  res1 (after source): " << (*res)("cell",0,9) << " " << (*res)("face",0,29) << std::endl;
#endif

  //print_vector(S_next_,res, "residual") ;
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void OverlandFlow::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if debug_flag
  std::cout << "Precon application:" << std::endl;
  std::cout << "  p0: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,3) << std::endl;
  std::cout << "  p1: " << (*u->data())("cell",0,9) << " " << (*u->data())("face",0,29) << std::endl;
#endif
  preconditioner_->ApplyInverse(*u->data(), Pu->data());
  //*Pu = *u ;
#if debug_flag
  std::cout << "  PC*p0: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,3) << std::endl;
  std::cout << "  PC*p1: " << (*Pu->data())("cell",0,9) << " " << (*Pu->data())("face",0,29) << std::endl;
#endif
};


// -----------------------------------------------------------------------------
// Compute a norm on (u,du)
// -----------------------------------------------------------------------------
double OverlandFlow::enorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) {
  double enorm_val = 0.0;
  Teuchos::RCP<const Epetra_MultiVector> pres_vec = u->data()->ViewComponent("cell", false);
  Teuchos::RCP<const Epetra_MultiVector> dpres_vec = du->data()->ViewComponent("cell", false);

  for (int lcv=0; lcv!=pres_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*dpres_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*pres_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
    //    printf("cell: %5i %14.7e %14.7e\n",lcv,(*(*dpres_vec)(0))[lcv],tmp);
  }
  

  Teuchos::RCP<const Epetra_MultiVector> fpres_vec = u->data()->ViewComponent("face", false);
  Teuchos::RCP<const Epetra_MultiVector> fdpres_vec = du->data()->ViewComponent("face", false);

  LINE(---) ;

  for (int lcv=0; lcv!=fpres_vec->MyLength(); ++lcv) {
    double tmp = abs((*(*fdpres_vec)(0))[lcv])/(atol_ + rtol_*abs((*(*fpres_vec)(0))[lcv]));
    enorm_val = std::max<double>(enorm_val, tmp);
    //    printf("face: %5i %14.7e %14.7e\n",lcv,(*(*fdpres_vec)(0))[lcv],tmp);
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  return enorm_val;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void OverlandFlow::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  std::cout << "Precon update at t = " << t << std::endl;
  // update state with the solution up.
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_);
  UpdateSecondaryVariables_(S_next_);

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_    ->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_);
  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetFieldData("upwind_overland_conductivity");

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  preconditioner_->CreateMFDstiffnessMatrices(*cond);
  preconditioner_->CreateMFDrhsVectors();

  // 2. Update local matrices diagonal with the accumulation terms.
  Teuchos::RCP<const CompositeVector> cell_volume =
    S_next_->GetFieldData("surface_cell_volume");

  Teuchos::RCP<const CompositeVector> pres = S_next_->GetFieldData("overland_pressure");
  Teuchos::RCP<const CompositeVector> rain = S_next_->GetFieldData("rainfall_rate");

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  int ncells = cell_volume->size("cell");
  for (int c=0; c!=ncells; ++c) {
    // accumulation term
    Acc_cells[c] += (*cell_volume)("cell",c) / h;
    Fc_cells[c] += (*pres)("cell",c) * (*cell_volume)("cell",c) / h;

    // source term
    //Fc_cells[c] += rhs_load_value(S_next_->time()) * (*cell_volume)("cell",0,c);
    Fc_cells[c] += (*rain)("cell",c) * (*cell_volume)("cell",c);
  }

  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();

  // currently the PC has p + z as the global RHS -- should correct
  // the global RHS to have the correct value (just p) prior to
  // calculating the Schur complement?  Maybe this is not necessary as
  // we're forming the Schur complement of the face-face system, so
  // the face RHS may not come into play for the operator (just the
  // Schur complement's RHS?)

  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);

  // Code to dump Schur complement to check condition number
  /*
    Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
    std::stringstream filename_s;
    filename_s << "schur_" << S_next_->cycle() << ".txt";
    //a  std::string filename = filename_s.str();
    EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
    std::cout << "updated precon " << S_next_->cycle() << std::endl;
  */

  preconditioner_->UpdateMLPreconditioner();
};

}  // namespace Flow
}  // namespace Amanzi
