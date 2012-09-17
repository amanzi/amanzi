/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "field_evaluator.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

// Richards is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Richards::fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  S_inter_->set_time(t_old);
  S_next_->set_time(t_new);
  double h = t_new - t_old;
  Teuchos::RCP<CompositeVector> u = u_new->data();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "----------------------------------------------------------------" << std::endl;
    *out_ << "Richards Residual calculation: T0 = " << t_old << " T1 = " << t_new << " H = " << h << std::endl;
    *out_ << "  p0: " << (*u)("cell",0,0) << " " << (*u)("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u)("cell",0,99) << " " << (*u)("face",0,497) << std::endl;
  }

  // pointer-copy temperature into state and update any auxilary data
  solution_to_state(u_new, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_, res);
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

    *out_ << "  res0 (after diffusion): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after diffusion): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;
  }

  // accumulation term
  AddAccumulation_(res);
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

    *out_ << "  res0 (after accumulation): " << (*res)("cell",0,0) << " " << (*res)("face",0,3) << std::endl;
    *out_ << "  res1 (after accumulation): " << (*res)("cell",0,99) << " " << (*res)("face",0,497) << std::endl;

    Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> darcy_flux = S_next_->GetFieldData("darcy_flux");

    *out_ << std::endl;
    //  *out_ << "  mass balance at 0: " << std::endl;

    //  *out_ << "      old water: " << (*wc0)("cell",0) << std::endl;
    //  *out_ << "      new water: " << (*wc1)("cell",0) << std::endl;

    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    wc0->mesh()->cell_get_faces_and_dirs(0, &faces, &dirs);
    double flux = 0.0;
    for (int lcv=0; lcv!=faces.size(); ++lcv) {
      //    *out_ << "      face " << faces[lcv] << ": " << dirs[lcv]*(*darcy_flux)("face",faces[lcv]) << std::endl;
      flux += dirs[lcv]*(*darcy_flux)("face",faces[lcv]);
    }
    //  *out_ << "    error: " << (*wc1)("cell",0) - (*wc0)("cell",0) + h*flux << std::endl;
    *out_ << "  mass balance error0: " << (*wc1)("cell",0) - (*wc0)("cell",0) + h*flux << std::endl;

    //  *out_ << std::endl;
    //  *out_ << "  mass balance at 99: " << std::endl;
    //  *out_ << "      old water: " << (*wc0)("cell",99) << std::endl;
    //  *out_ << "      new water: " << (*wc1)("cell",99) << std::endl;

    faces.clear();
    dirs.clear();
    wc0->mesh()->cell_get_faces_and_dirs(99, &faces, &dirs);
    flux = 0.0;
    for (int lcv=0; lcv!=faces.size(); ++lcv) {
      //    *out_ << "      face " << faces[lcv] << ": " << dirs[lcv]*(*darcy_flux)("face",faces[lcv]) << std::endl;
      flux += dirs[lcv]*(*darcy_flux)("face",faces[lcv]);
    }
    //  *out_ << "    error: " << (*wc1)("cell",99) - (*wc0)("cell",99) + h*flux << std::endl;
    *out_ << "  mass balance error1: " << (*wc1)("cell",99) - (*wc0)("cell",99) + h*flux << std::endl;
  }

};

// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
void Richards::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Precon application:" << std::endl;
    *out_ << "  p0: " << (*u->data())("cell",0,0) << " " << (*u->data())("face",0,3) << std::endl;
    *out_ << "  p1: " << (*u->data())("cell",0,99) << " " << (*u->data())("face",0,497) << std::endl;
  }

  Teuchos::RCP<const CompositeVector> pres = u->data();
  for (int c=0; c!=pres->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*pres)("cell",c))) {
      *out_ << "Cutting time step due to NaN in cell residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=pres->size("face",false); ++f) {
    if (boost::math::isnan<double>((*pres)("face",f))) {
      *out_ << "Cutting time step due to NaN in face residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

  preconditioner_->ApplyInverse(*u->data(), Pu->data());

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

  *out_ << "  PC*p0: " << (*Pu->data())("cell",0,0) << " " << (*Pu->data())("face",0,3) << std::endl;
  *out_ << "  PC*p1: " << (*Pu->data())("cell",0,99) << " " << (*Pu->data())("face",0,497) << std::endl;
  }

  Teuchos::RCP<const CompositeVector> ppres = Pu->data();
  for (int c=0; c!=ppres->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*ppres)("cell",c))) {
      // dump the schur complement
      Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
      std::stringstream filename_s;
      filename_s << "schur_" << S_next_->cycle() << ".txt";
      EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
      *out_ << "updated precon " << S_next_->cycle() << std::endl;

      // print the rel perm
      Teuchos::RCP<const CompositeVector> num_rel_perm = S_next_->GetFieldData("numerical_rel_perm");
      Teuchos::RCP<const CompositeVector> rel_perm = S_next_->GetFieldData("relative_permeability");
      *out_ << "REL PERM: " << std::endl;
      rel_perm->Print(*out_);
      *out_ << std::endl;
      *out_ << "UPWINDED REL PERM: " << std::endl;
      num_rel_perm->Print(*out_);


      // throw
      *out_ << "Cutting time step due to NaN in PC'd cell residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=ppres->size("face",false); ++f) {
    if (boost::math::isnan<double>((*ppres)("face",f))) {
      *out_ << "Cutting time step due to NaN in PC'd face residual." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }

};


// -----------------------------------------------------------------------------
// Compute a norm on (u,du)
// -----------------------------------------------------------------------------
  /*
double Richards::enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  Teuchos::RCP<const CompositeVector> pres = u->data();
  Teuchos::RCP<const CompositeVector> dpres = du->data();

  // error in mass conservation, which is what we really care about.  This
  // would need serious work to be optimized, this is crap having to duplicate
  // everything.
  double enorm_mass = 0.0;
  PK::solution_to_state(u, S_next_);
  double h = S_next_->time() - S_inter_->time();

  bool update = S_next_->GetFieldEvaluator("relative_permeability")->HasFieldChanged(S_next_.ptr(), name_);
  if (update) UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("numerical_rel_perm", name_);

  matrix_->CreateMFDstiffnessMatrices(*rel_perm);
  Teuchos::RCP<CompositeVector> darcy_flux =
    S_next_->GetFieldData("darcy_flux", name_);
  matrix_->DeriveFlux(*pres, darcy_flux);
  AddGravityFluxesToVector_(S_next_, darcy_flux);
  S_next_->GetFieldEvaluator("water_content")->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");

  for (int c=0; c!=wc1->size("cell"); ++c) {
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    wc0->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    double flux = 0.0;
    for (int lcv=0; lcv!=faces.size(); ++lcv) {
      flux += dirs[lcv]*(*darcy_flux)("face",faces[lcv]);
    }
    double tmp = (*wc1)("cell",c) - (*wc0)("cell",c) + h*flux;
    tmp = tmp / (mass_atol_ + mass_rtol_ * (*wc0)("cell",c));
    enorm_mass = std::max<double>(enorm_mass, tmp);
  }


  *out_ << "ENORM (cell, face, mass): " << std::scientific << enorm_val_cell
            << " / " << std::scientific << enorm_val_face
            << " / " << std::scientific << enorm_mass << std::endl;

    double enorm_val = std::max<double>(enorm_val_cell, enorm_val_face);
  //  enorm_val = std::max<double>(enorm_val, enorm_mass);

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  return enorm_val;
};
  */

// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Richards::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

  *out_ << "Precon update at t = " << t << std::endl;
  }

  // update state with the solution up.
  S_next_->set_time(t);
  PKDefaultBase::solution_to_state(up, S_next_);

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_next_->GetFieldData("numerical_rel_perm");

  preconditioner_->CreateMFDstiffnessMatrices(*rel_perm);
  preconditioner_->CreateMFDrhsVectors();
  AddGravityFluxes_(S_next_, preconditioner_);


  // update with accumulation terms
  // -- update the accumulation derivatives
  S_next_->GetFieldEvaluator("water_content")
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);

  // -- get the accumulation deriv
  Teuchos::RCP<const CompositeVector> dwc_dp =
      S_next_->GetFieldData("dwater_content_dpressure");
  Teuchos::RCP<const CompositeVector> pres =
      S_next_->GetFieldData(key_);

  /*
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {

  // update with accumulation terms
  Teuchos::RCP<const double> p_atm = S_next_->GetScalarData("atmospheric_pressure");
  Teuchos::RCP<const CompositeVector> temp = S_next_->GetFieldData("temperature");
  Teuchos::RCP<const CompositeVector> poro = S_next_->GetFieldData("porosity");

  Teuchos::RCP<const CompositeVector> mol_frac_gas =
    S_next_->GetFieldData("mol_frac_gas");
  Teuchos::RCP<const CompositeVector> n_gas = S_next_->GetFieldData("molar_density_gas");
  Teuchos::RCP<const CompositeVector> sat_gas = S_next_->GetFieldData("saturation_gas");

  Teuchos::RCP<const CompositeVector> n_liq = S_next_->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> sat_liq = S_next_->GetFieldData("saturation_liquid");

  Teuchos::RCP<const CompositeVector> n_ice = S_next_->GetFieldData("molar_density_ice");
  Teuchos::RCP<const CompositeVector> sat_ice = S_next_->GetFieldData("saturation_ice");

  Teuchos::RCP<const CompositeVector> cell_volume = S_next_->GetFieldData("cell_volume");

  int c = 0;
  double p = (*pres)("cell",c);
  double T = (*temp)("cell",c);
  double phi = (*poro)("cell",c);

  *out_ << "    p =" << p << std::endl;
  *out_ << "    T =" << T << std::endl;
  *out_ << "    phi =" << phi << std::endl;
  *out_ << "    cv =" << (*cell_volume)("cell",c) << std::endl;
  *out_ << "   res3 (0) =" << (*dwc_dp)("cell",c) / phi / (*cell_volume)("cell",c) << std::endl;

  c = 99;
  p = (*pres)("cell",c);
  T = (*temp)("cell",c);
  phi = (*poro)("cell",c);
  *out_ << "    p =" << p << std::endl;
  *out_ << "    T =" << T << std::endl;
  *out_ << "    phi =" << phi << std::endl;
  *out_ << "    cv =" << (*cell_volume)("cell",c) << std::endl;
  *out_ << "   res3 (99) =" << (*dwc_dp)("cell",c) / phi / (*cell_volume)("cell",c) << std::endl;
  }
  */

  // -- get the matrices/rhs that need updating
  std::vector<double>& Acc_cells = preconditioner_->Acc_cells();
  std::vector<double>& Fc_cells = preconditioner_->Fc_cells();

  int ncells = pres->size("cell");
  for (int c=0; c!=ncells; ++c) {
    Acc_cells[c] += (*dwc_dp)("cell",c) / h;
    //    Fc_cells[c] += (*dwc_dp)("cell",c) / h * (*pres)("cell",c);
  }

  preconditioner_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  preconditioner_->AssembleGlobalMatrices();
  preconditioner_->ComputeSchurComplement(bc_markers_, bc_values_);

  // Code to dump Schur complement to check condition number
  /*
  Teuchos::RCP<Epetra_FECrsMatrix> sc = preconditioner_->Schur();
  std::stringstream filename_s;
  filename_s << "schur_" << S_next_->cycle() << ".txt";
  //a  std::string filename = filename_s.str();
  EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  *out_ << "updated precon " << S_next_->cycle() << std::endl;
  */

  preconditioner_->UpdatePreconditioner();
};

}  // namespace Flow
}  // namespace Amanzi



