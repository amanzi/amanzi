/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>
#include <boost/math/tools/roots.hpp>

#include "OperatorDefs.hh"
#include "OperatorDiffusionFVwithGravity.hh"
#include "Op.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Op_Face_Cell.hh"
#include "Operator_Cell.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization completes initialization of the base class.
****************************************************************** */
void OperatorDiffusionFVwithGravity::InitDiffusion_(Teuchos::ParameterList& plist)
{
  gravity_term_ = Teuchos::rcp(new CompositeVector(*transmissibility_));
  g_.set(mesh_->space_dimension(), 0.0);
}


/* ******************************************************************
* Setup methods: scalar coefficients
****************************************************************** */
void OperatorDiffusionFVwithGravity::Setup(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K)
{
  K_ = K;
  
  if (!transmissibility_initialized_) {
    ComputeTransmissibility_(&g_, gravity_term_);
  }
}


/* ******************************************************************
* Setup methods: krel and deriv -- must be called after calling a
* setup with K absolute
****************************************************************** */
void OperatorDiffusionFVwithGravity::Setup(const Teuchos::RCP<const CompositeVector>& k,
                                           const Teuchos::RCP<const CompositeVector>& dkdp)
{
  k_ = k;
  dkdp_ = dkdp;

  if (k_ != Teuchos::null) {
    ASSERT(k_->HasComponent("face"));
    // NOTE: it seems that Amanzi passes in a cell based kr which is then
    // ignored, and assumed = 1.  This seems dangerous to me. --etc
    // ASSERT(!k_->HasComponent("cell"));
  }
  if (dkdp_ != Teuchos::null) {
    ASSERT(dkdp_->HasComponent("cell"));
    ASSERT(dkdp_->HasComponent("face"));
  }

  // verify that mass matrices were initialized.
  if (!transmissibility_initialized_) {
    ComputeTransmissibility_(&g_, gravity_term_);
  }
}


/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void OperatorDiffusionFVwithGravity::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                                    const Teuchos::Ptr<const CompositeVector>& u)
{
  OperatorDiffusionFV::UpdateMatrices(flux, u);

  // populating right-hand side
  if (!exclude_primary_terms_) {
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;

    const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
    const Epetra_MultiVector& gravity_face = *gravity_term_->ViewComponent("face", true);
    Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

    // preparing upwind data
    Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
    if (k_ != Teuchos::null) {
      if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
    }

    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        if (bc_model[f] == OPERATOR_BC_NEUMANN) continue;
        rhs_cell[0][c] -= dirs[n] * gravity_face[0][f] * (k_face.get() ? (*k_face)[0][f] : 1.0);
      }
    }
  }
}


/* ******************************************************************
* Special implementation of boundary conditions.
**********************************;******************************** */
void OperatorDiffusionFVwithGravity::ApplyBCs(bool primary, bool eliminate)
{
  OperatorDiffusionFV::ApplyBCs(primary, eliminate);

  Epetra_MultiVector* gravity_face = &*gravity_term_->ViewComponent("face", true);
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == OPERATOR_BC_NEUMANN) {
      (*gravity_face)[0][f] = 0.0;
    }
  }
}


/* ******************************************************************
* Calculate flux from cell-centered data
****************************************************************** */
void OperatorDiffusionFVwithGravity::UpdateFlux(
    const CompositeVector& solution, CompositeVector& darcy_mass_flux)
{
  OperatorDiffusionFV::UpdateFlux(solution, darcy_mass_flux);

  // add gravity term
  Epetra_MultiVector& flux = *darcy_mass_flux.ViewComponent("face", false);
  const Teuchos::Ptr<const Epetra_MultiVector> Krel_face =
      k_.get() ? k_->ViewComponent("face", false).ptr() : Teuchos::null;

  if (Krel_face.get()) {
    flux.Multiply(1.0, *Krel_face, *gravity_term_->ViewComponent("face", false), 1.0);
  } else {
    flux.Update(1.0, *gravity_term_->ViewComponent("face", false), 1.0);
  }
}


/* ******************************************************************
* Computation of a local submatrix of the analytical Jacobian 
* (its nonlinear part) on face f.
****************************************************************** */
void OperatorDiffusionFVwithGravity::ComputeJacobianLocal_(
    int mcells, int f, int face_dir, int Krel_method,
    int bc_model_f, double bc_value_f,
    double *pres, double *dkdp_cell,
    WhetStone::DenseMatrix& Jpp)
{
  const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
  const Teuchos::Ptr<Epetra_MultiVector> gravity_face =
      gravity_term_.get() ? gravity_term_->ViewComponent("face", true).ptr() : Teuchos::null;

  double dKrel_dp[2];
  double dpres;

  if (mcells == 2) {
    dpres = pres[0] - pres[1];  // + grn;
    if (Krel_method == OPERATOR_LITTLE_K_UPWIND) {
      double flux0to1;
      flux0to1 = trans_face[0][f] * dpres;
      if (gravity_face.get()) flux0to1 += face_dir * (*gravity_face)[0][f];
      if (flux0to1  > OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dkdp_cell[0];
        dKrel_dp[1] = 0.0;
      } else if (flux0to1 < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dkdp_cell[1];
      } else if (fabs(flux0to1) < OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.5 * dkdp_cell[0];
        dKrel_dp[1] = 0.5 * dkdp_cell[1];
      }
    } else if (Krel_method == OPERATOR_UPWIND_ARITHMETIC_AVERAGE) {
      dKrel_dp[0] = 0.5 * dkdp_cell[0];
      dKrel_dp[1] = 0.5 * dkdp_cell[1];
    } else {
      ASSERT(0);
    }

    Jpp(0, 0) = trans_face[0][f] * dpres * dKrel_dp[0];
    Jpp(0, 1) = trans_face[0][f] * dpres * dKrel_dp[1];

    if (gravity_face.get()) {
      Jpp(0,0) += face_dir * (*gravity_face)[0][f] * dKrel_dp[0];
      Jpp(0,1) += face_dir * (*gravity_face)[0][f] * dKrel_dp[1];
    }
    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_model_f == OPERATOR_BC_DIRICHLET) {                   
      pres[1] = bc_value_f;
      dpres = pres[0] - pres[1];  // + grn;
      Jpp(0, 0) = trans_face[0][f] * dpres * dkdp_cell[0];
      if (gravity_face.get())
          Jpp(0,0) += face_dir * (*gravity_face)[0][f] * dkdp_cell[0];
    } else {
      Jpp(0, 0) = 0.0;
    }
  }
}


/* ******************************************************************
* Return value of the gravity flux on the given face f.
****************************************************************** */
double OperatorDiffusionFVwithGravity::ComputeGravityFlux(int f) const
{
  const Epetra_MultiVector& gravity_face = *gravity_term_->ViewComponent("face", true); 
  return gravity_face[0][f];
}

}  // namespace Operators
}  // namespace Amanzi


