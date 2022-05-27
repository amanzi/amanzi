/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Teuchos_RCP.hpp"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "InverseFactory.hh"
#include "PDE_Accumulation.hh"
#include "PDE_DiffusionFactory.hh"
#include "MFD3D_Diffusion.hh"
#include "nlfv.hh"
#include "Tensor.hh"

#include "TransportDefs.hh"
#include "Transport_PK.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
* Calculate dispersive tensor from given Darcy fluxes. The flux is
* assumed to be scaled by face area.
******************************************************************* */
void Transport_PK::CalculateDispersionTensor_(
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  if (!flag_dispersion_) {
    D_.clear();
    return;
  }

  D_.resize(ncells_owned);
  for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);

  const auto& flowrate = *S_->Get<CompositeVector>(vol_flowrate_key_).ViewComponent("face", true);

  AmanziGeometry::Point velocity(dim);
  WhetStone::MFD3D_Diffusion mfd3d(mesh_);
  WhetStone::Polynomial poly(dim, 1);

  for (int c = 0; c < ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    std::vector<WhetStone::Polynomial> flux(nfaces);
    for (int n = 0; n < nfaces; n++) {
      flux[n].Reshape(dim, 0);
      flux[n](0) = flowrate[0][faces[n]];
    }
    mfd3d.L2Cell(c, flux, flux, NULL, poly);

    for (int k = 0; k < dim; ++k) velocity[k] = poly(k + 1);
    D_[c] = mdm_->second[(*mdm_->first)[c]]->mech_dispersion(
        velocity, axi_symmetry_[c], saturation[0][c], porosity[0][c]);
  }
}


/* *******************************************************************
* Calculate diffusion tensor and add it to the dispersion tensor.
******************************************************************* */
void Transport_PK::CalculateDiffusionTensor_(
    double md, int phase, 
    const Epetra_MultiVector& porosity, const Epetra_MultiVector& saturation)
{
  if (D_.size() == 0) {
    D_.resize(ncells_owned);
    for (int c = 0; c < ncells_owned; c++) D_[c].Init(dim, 1);
  }

  for (int mb = 0; mb < mat_properties_.size(); mb++) {
    Teuchos::RCP<MaterialProperties> spec = mat_properties_[mb]; 

    std::vector<AmanziMesh::Entity_ID> block;
    for (int r = 0; r < (spec->regions).size(); r++) {
      std::string region = (spec->regions)[r];
      mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

      AmanziMesh::Entity_ID_List::iterator c;
      if (phase == TRANSPORT_PHASE_LIQUID) {
        for (c = block.begin(); c != block.end(); c++) {
          D_[*c] += md * spec->tau[phase] * porosity[0][*c] * saturation[0][*c];
        }
      } else if (phase == TRANSPORT_PHASE_GAS) {
        for (c = block.begin(); c != block.end(); c++) {
          D_[*c] += md * spec->tau[phase] * porosity[0][*c] * (1.0 - saturation[0][*c]);
        }
      }
    }
  }
}


/* ******************************************************************
* Check all phases for the given name.
****************************************************************** */
int Transport_PK::FindDiffusionValue(const std::string& tcc_name, double* md, int* phase)
{
  for (int i = 0; i < TRANSPORT_NUMBER_PHASES; i++) {
    if (diffusion_phase_[i] == Teuchos::null) continue;
    int ok = diffusion_phase_[i]->FindDiffusionValue(tcc_name, md);
    if (ok == 0) {
      *phase = i;
      return 0;
    }
  }

  *md = 0.0;
  *phase = -1;
  return -1;
}


/* ******************************************************************
* Find direction of axi-symmetry for Lichtner-Kelkar-Robinson model
****************************************************************** */
void Transport_PK::CalculateAxiSymmetryDirection()
{
  axi_symmetry_.resize(ncells_owned, -1);
  if (S_->HasRecord(permeability_key_) && dim == 3) {
    const auto& perm = *S_->Get<CompositeVector>(permeability_key_).ViewComponent("cell");

    for (int c = 0; c < ncells_owned; ++c) {
      int k = -1;
      if (perm[0][c] != perm[1][c] && perm[1][c] == perm[2][c]) {
        k = 0;
      } else if (perm[1][c] != perm[2][c] && perm[2][c] == perm[0][c]) {
        k = 1;
      } else if (perm[2][c] != perm[0][c] && perm[0][c] == perm[1][c]) {
        k = 2;
      } 
      axi_symmetry_[c] = k;
    }
  }
}


/* ******************************************************************
* Two-phase solver
****************************************************************** */
Teuchos::RCP<Operators::Operator> Transport_PK::DispersionSolver(
    const Epetra_MultiVector& tcc_prev,
    const Epetra_MultiVector& tcc_next,
    double t_old, double t_new, int comp0)
{
  CalculateDispersionTensor_(*transport_phi, *ws);

  int i0 = (comp0 >= 0) ? comp0 : 0;

  int num_components = tcc_prev.NumVectors();
  double dt_MPC = t_new - t_old;

  auto bc_dummy = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // create the Dispersion and Accumulation operators
  Teuchos::ParameterList& op_list = tp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");
  Operators::PDE_DiffusionFactory opfactory;
  Teuchos::RCP<Operators::PDE_Diffusion> op1 = opfactory.Create(op_list, mesh_, bc_dummy);

  op1->SetBCs(bc_dummy, bc_dummy);
  auto op = op1->global_operator();
  auto op2 = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op));

  // Create the preconditioner and solver.
  //
  // This implementation re-create the Ops each timestep.  This could be
  // updated to re-use them, but would require mass matrices to be
  // recalculated. --etc
  auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
      dispersion_preconditioner, *preconditioner_list_,
      dispersion_solver, *linear_solver_list_, true);
  inv_list.setName(dispersion_preconditioner);
  op->set_inverse_parameters(inv_list);
  op->InitializeInverse();

  const CompositeVectorSpace& cvs = op1->global_operator()->DomainMap();
  CompositeVector sol(cvs), factor(cvs), factor0(cvs), source(cvs);

  int phase, num_itrs(0);
  double md_change, md_old(0.0), md_new, residual(0.0);

  // Disperse and diffuse aqueous components
  // -- if we reach this point, we instantiate diffusion operator 
  //    regadless of flags, since MPC may need a properly formed operator 
  //    (e.g. with face DOFs)
  for (int i = i0; i < num_aqueous; i++) {
    FindDiffusionValue(component_names_[i], &md_new, &phase);
    md_change = md_new - md_old;
    md_old = md_new;

    if (md_change != 0.0 || D_.size() == 0) {
      CalculateDiffusionTensor_(md_change, phase, *transport_phi, *ws);
    }

    // set the initial guess
    Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) {
      sol_cell[0][c] = tcc_next[i][c];
    }
    if (sol.HasComponent("face")) {
      sol.ViewComponent("face")->PutScalar(0.0);
    }

    op->Init();
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Dptr = Teuchos::rcpFromRef(D_);
    op1->Setup(Dptr, Teuchos::null, Teuchos::null);
    op1->UpdateMatrices(Teuchos::null, Teuchos::null);

    // boundary conditions
    std::vector<int>& bc_model = bc_dummy->bc_model();
    std::vector<double>& bc_value = bc_dummy->bc_value();
    ComputeBCs_(bc_model, bc_value, i);

    // add accumulation term
    Epetra_MultiVector& fac = *factor.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) {
      fac[0][c] = (*phi)[0][c] * (*ws)[0][c];
    }
    op2->AddAccumulationDelta(sol, factor, factor, dt_MPC, "cell");
    op1->ApplyBCs(true, true, true);
    if (comp0 >= 0) return op;

    op->ComputeInverse();
    CompositeVector& rhs = *op->rhs();
    int ierr = op->ApplyInverse(rhs, sol);

    if (ierr < 0) {
      Errors::Message msg("TransportExplicit_PK solver failed with message: \"");
      msg << op->returned_code_string() << "\"";
      Exceptions::amanzi_throw(msg);
    }

    residual += op->residual();
    num_itrs += op->num_itrs();

    for (int c = 0; c < ncells_owned; c++) {
      tcc_next[i][c] = sol_cell[0][c];
    }
  }

  // Diffuse gaseous components. We ignore dispersion 
  // tensor (D is reset). Inactive cells (sat_l[c] = 1 and D_[c] = 0) 
  // are treated with a hack of the accumulation term.
  D_.clear();
  md_old = 0.0;
  for (int i = std::max(i0, num_aqueous); i < num_components; i++) {
    FindDiffusionValue(component_names_[i], &md_new, &phase);
    md_change = md_new - md_old;
    md_old = md_new;

    if (md_change != 0.0 || i == num_aqueous) {
      CalculateDiffusionTensor_(md_change, phase, *transport_phi, *ws);
    }

    // set initial guess
    Epetra_MultiVector& sol_cell = *sol.ViewComponent("cell");
    for (int c = 0; c < ncells_owned; c++) {
      sol_cell[0][c] = tcc_next[i][c];
    }
    if (sol.HasComponent("face")) {
      sol.ViewComponent("face")->PutScalar(0.0);
    }

    op->Init();
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Dptr = Teuchos::rcpFromRef(D_);
    op1->Setup(Dptr, Teuchos::null, Teuchos::null);
    op1->UpdateMatrices(Teuchos::null, Teuchos::null);

    // add boundary conditions and sources for gaseous components
    std::vector<int>& bc_model = bc_dummy->bc_model();
    std::vector<double>& bc_value = bc_dummy->bc_value();
    ComputeBCs_(bc_model, bc_value, i);

    Epetra_MultiVector& rhs_cell = *op->rhs()->ViewComponent("cell");
    ComputeSources_(t_new, 1.0, rhs_cell, tcc_prev, i, i);
    op1->ApplyBCs(true, true, true);

    // add accumulation term
    Epetra_MultiVector& fac1 = *factor.ViewComponent("cell");
    Epetra_MultiVector& fac0 = *factor0.ViewComponent("cell");

    for (int c = 0; c < ncells_owned; c++) {
      fac1[0][c] = (*phi)[0][c] * (1.0 - (*ws)[0][c]);
      fac0[0][c] = (*phi)[0][c] * (1.0 - (*ws_prev)[0][c]);
      if ((*ws)[0][c] == 1.0) fac1[0][c] = 1.0;  // hack so far
    }
    op2->AddAccumulationDelta(sol, factor0, factor, dt_MPC, "cell");
    if (comp0 >= 0) return op;

    op->ComputeInverse();
    CompositeVector& rhs = *op->rhs();
    int ierr = op->ApplyInverse(rhs, sol);

    if (ierr < 0) {
      Errors::Message msg("TransportExplicit_PK solver failed with message: \"");
      msg << op->returned_code_string() << "\"";
      Exceptions::amanzi_throw(msg);
    }

    residual += op->residual();
    num_itrs += op->num_itrs();

    for (int c = 0; c < ncells_owned; c++) {
      tcc_next[i][c] = sol_cell[0][c];
    }
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "dispersion solver (" << dispersion_solver
               << ") ||r||=" << residual / num_components
               << " itrs=" << num_itrs / num_components << std::endl;
  }

  return Teuchos::null;
}

}  // namespace Transport
}  // namespace Amanzi



