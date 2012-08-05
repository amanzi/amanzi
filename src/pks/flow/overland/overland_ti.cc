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

#if 0
}}
#endif

#define DEBUG_FLAG 0
#define DEBUG_RES_FLAG 0
#define UPDATE_FOR_REAL 0
#define DEBUG_ERROR_FLAG 0

// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void OverlandFlow::fun( double t_old,
                        double t_new,
                        Teuchos::RCP<TreeVector> u_old,
                        Teuchos::RCP<TreeVector> u_new,
                        Teuchos::RCP<TreeVector> g ) {
  niter_++;

  S_inter_->set_time(t_old);
  S_next_ ->set_time(t_new);

  Teuchos::RCP<CompositeVector> u = u_new->data();
#if DEBUG_FLAG
  std::cout << "OverlandFlow Residual calculation:" << std::endl;
  std::cout << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,0) << std::endl;
  std::cout << "  p1: " << (*u)("cell",0,9) << " " << (*u)("face",0,29) << std::endl;
#endif
  //print_vector(S_next_,u,"fun") ;

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // ensure postive pressure
  // double minval = 0.0;
  // S_next_->GetFieldData("overland_pressure")->ViewComponent("cell", false)->MinValue(&minval);
  // if (minval < 0.0) {
  //   // cannot handle negative pressures, punt
  //   std::cout << "Cutting time step due negative pressure in residual function." << std::endl;
  //   Errors::Message m("Cut time step");
  //   Exceptions::amanzi_throw(m);
  // }

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

  // update the preconditioner
  UpdateBoundaryConditionsNoElev_(S_next_);
  update_precon_for_real(t_new, u_new, t_new - t_old);

#if DEBUG_RES_FLAG
  std::cout << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;
  std::cout << "  res1 (after diffusion): " << (*res)("cell",0,9) << " " << (*res)("face",0,29) << std::endl;
#endif

  // accumulation term
  AddAccumulation_(res);
#if DEBUG_RES_FLAG
  std::cout << "  res0 (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;
  std::cout << "  res1 (after accumulation): " << (*res)("cell",0,9) << " " << (*res)("face",0,29) << std::endl;
#endif

  // add rhs load value
  AddLoadValue_(S_next_,res);
#if DEBUG_RES_FLAG
  std::cout << "  res0 (after source): " << (*res)("cell",0,0) << " " << (*res)("face",0,0) << std::endl;
  std::cout << "  res1 (after source): " << (*res)("cell",0,9) << " " << (*res)("face",0,29) << std::endl;
#endif

  //print_vector(S_next_,res, "residual") ;
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void OverlandFlow::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {

#if DEBUG_FLAG
  std::cout << "Precon application:" << std::endl;
  std::cout << "  p0: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,0) << std::endl;
  std::cout << "  p1: " << (*u->data())("cell",0,9) << " " << (*u->data())("face",0,29) << std::endl;
#endif

  // check for nans in residual
  Teuchos::RCP<const CompositeVector> res = u->data();
  for (int c=0; c!=res->size("cell"); ++c) {
    if (boost::math::isnan<double>((*res)("cell",c))) {
      int mypid = S_next_->Mesh("surface")->get_comm()->MyPID();
      std::cout << "Cutting time step due to NaN in cell residual on proc " << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=res->size("face"); ++f) {
    if (boost::math::isnan<double>((*res)("face",f))) {
      int mypid = S_next_->Mesh("surface")->get_comm()->MyPID();
      std::cout << "Cutting time step due to NaN in face residual on proc " << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  // check for nans in preconditioned residual
  Teuchos::RCP<const CompositeVector> Pres = Pu->data();
  for (int c=0; c!=Pres->size("cell"); ++c) {
    if (boost::math::isnan<double>((*Pres)("cell",c))) {
      int mypid = S_next_->Mesh("surface")->get_comm()->MyPID();
      std::cout << "Cutting time step due to NaN in PC'd cell residual on proc " << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=Pres->size("face"); ++f) {
    if (boost::math::isnan<double>((*Pres)("face",f))) {
      int mypid = S_next_->Mesh("surface")->get_comm()->MyPID();
      std::cout << "Cutting time step due to NaN in PC'd face residual on proc " << mypid << "." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

  //*Pu = *u ;
#if DEBUG_FLAG
  std::cout << "  PC*p0: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,0) << std::endl;
  std::cout << "  PC*p1: " << (*Pu->data())("cell",0,9) << " " << (*Pu->data())("face",0,29) << std::endl;
#endif
};


// -----------------------------------------------------------------------------
// Compute a norm on (u,du)
// -----------------------------------------------------------------------------
double OverlandFlow::enorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) {
  Teuchos::RCP<const CompositeVector> pres = u->data();
  Teuchos::RCP<const CompositeVector> dpres = du->data();

  Teuchos::RCP<const CompositeVector> cell_volume =
    S_next_->GetFieldData("surface_cell_volume");
  Teuchos::RCP<const CompositeVector> slope =
    S_next_->GetFieldData("slope_magnitude");
  Teuchos::RCP<const CompositeVector> relperm =
    S_next_->GetFieldData("upwind_overland_conductivity");


#if DEBUG_FLAG
  std::cout << "enorm:" << std::endl;
  std::cout << "  c0 (p, dp): " << (*pres)("cell",0,0) << " " << (*dpres)("cell",0,0) << std::endl;
  std::cout << "  f0 (p, dp): " << (*pres)("face",0,0) << " " << (*dpres)("face",0,0) << std::endl;
  std::cout << "  c9 (p, dp): " << (*pres)("cell",0,9) << " " << (*dpres)("cell",0,9) << std::endl;
  std::cout << " f29 (p, dp): " << (*pres)("face",0,29) << " " << (*dpres)("face",0,29) << std::endl;
#endif

  double enorm_val_cell = 0.0;
  for (int c=0; c!=pres->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*dpres)("cell",c))) {
      std::cout << "Cutting time step due to NaN in correction." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }

    double tmp = abs((*dpres)("cell",c)) / (atol_ + rtol_ * abs((*pres)("cell",c)));
    enorm_val_cell = std::max<double>(enorm_val_cell, tmp);
    //    printf("cell: %5i %14.7e %14.7e\n",lcv,(*(*dpres_vec)(0))[lcv],tmp);
  }

  double enorm_val_face = 0.0;
  for (int f=0; f!=pres->size("face",false); ++f) {
    if (boost::math::isnan<double>((*dpres)("face",f))) {
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }

    double tmp = abs((*dpres)("face",f)) / (atol_ + rtol_ * abs((*pres)("face",f)));
    enorm_val_face = std::max<double>(enorm_val_face, tmp);
    //    printf("face: %5i %14.7e %14.7e\n",lcv,(*(*fdpres_vec)(0))[lcv],tmp);
  }


  //  std::cout.precision(15);
  //  std::cout << "enorm val (cell, face): " << std::scientific << enorm_val_cell
  //            << " / " << std::scientific << enorm_val_face << std::endl;

  double enorm_val = std::max<double>(enorm_val_cell, enorm_val_face);
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
  // fooled you.
}


void OverlandFlow::update_precon_for_real(double t, Teuchos::RCP<const TreeVector> up, double h) {
#if DEBUG_FLAG
  std::cout << "Precon update at t = " << t << std::endl;
  std::cout << "  p0: " << (*up->data())("cell",0,0) << " " << (*up->data())("face",0,0) << std::endl;
  std::cout << "  p1: " << (*up->data())("cell",0,9) << " " << (*up->data())("face",0,29) << std::endl;
#endif


  // update state with the solution up.
#ifdef UPDATE_FOR_REAL
  S_next_->set_time(t);
  PK::solution_to_state(up, S_next_);
  UpdateSecondaryVariables_(S_next_);
#endif

#ifdef UPDATE_FOR_REAL
  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_    ->Compute(S_next_->time());
  UpdateBoundaryConditionsNoElev_(S_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_);
#endif

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
  Teuchos::RCP<const CompositeVector> slope_mag = S_next_->GetFieldData("slope_magnitude");

  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();
  int ncells = cell_volume->size("cell");
  for (int c=0; c!=ncells; ++c) {
    // accumulation term
    Acc_cells[c] += (*cell_volume)("cell",c) / h;
    Fc_cells[c] += (*pres)("cell",c) * (*cell_volume)("cell",c) / h;

    // source term
    //Fc_cells[c] += rhs_load_value(S_next_->time()) * (*cell_volume)("cell",0,c);
    //Fc_cells[c] += (*rain)("cell",c) * (*cell_volume)("cell",c);
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
    // Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
    // std::stringstream filename_s;
    // filename_s << "schur_" << S_next_->cycle() << ".txt";
    // //a  std::string filename = filename_s.str();
    // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
    // std::cout << "updated precon " << S_next_->cycle() << std::endl;

  preconditioner_->UpdateMLPreconditioner();

  //  test_precon(t,up,h);
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
      AmanziGeometry::Point point = S_next_->Mesh("surface")->cell_centroid(c);
      if (error > 1e-5) {
        std::cout << "Bad error at cell: " << c << std::endl;
        AmanziMesh::Entity_ID_List faces;
        S_next_->Mesh("surface")->cell_get_faces(c, &faces);
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
