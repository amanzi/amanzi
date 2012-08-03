/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Neil Carlson (version 1)
  Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */

#include "field_evaluator.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Diffusion term, div K grad p
// -------------------------------------------------------------
void Richards::ApplyDiffusion_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& g) {

  // update the rel perm on cells.
  S->GetFieldEvaluator("relative_permeability")->HasFieldChanged(S.ptr(), "richards_pk");

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> rel_perm =
    S->GetFieldData("numerical_rel_perm", "flow");

  // update the stiffness matrix
  matrix_->CreateMFDstiffnessMatrices(*rel_perm);
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_(S, matrix_);
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  matrix_->ComputeNegativeResidual(*pres, g);
};


// -------------------------------------------------------------
// Accumulation of water term du/dt
// -------------------------------------------------------------
void Richards::AddAccumulation_(const Teuchos::RCP<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the water content at both the old and new times.
  S_next_->GetFieldEvaluator("water_content")->HasFieldChanged(S_next_.ptr(), "richards_pk");
  S_inter_->GetFieldEvaluator("water_content")->HasFieldChanged(S_inter_.ptr(), "richards_pk");

  // get these fields
  Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");

  // Water content only has cells, while the residual has cells and faces, so
  // this requires a little care.
  g->ViewComponent("cell",false)->Update(1.0/dt, *wc1->ViewComponent("cell",false),
          -1.0/dt, *wc0->ViewComponent("cell",false), 1.0);
};


// -------------------------------------------------------------
// Convert abs perm vector to tensor.
// -------------------------------------------------------------
void Richards::SetAbsolutePermeabilityTensor_(const Teuchos::RCP<State>& S) {
  // currently assumes isotropic perm, should be updated
  S->GetFieldEvaluator("permeability")->HasFieldChanged(S.ptr(), "richards_pk");
  Teuchos::RCP<const CompositeVector> perm = S->GetFieldData("permeability");
  int ncells = perm->size("cell");
  int ndofs = perm->num_dofs("cell");

  if (ndofs == 1) { // isotropic
    for (int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = (*perm)("cell",c);
    }
  } else if (ndofs == 2 && S->Mesh()->space_dimension() == 3) {
    // horizontal and vertical perms
    for (int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = (*perm)("cell",0,c);
      (*K_)[c](1, 1) = (*perm)("cell",0,c);
      (*K_)[c](2, 2) = (*perm)("cell",1,c);
    }
  } else if (ndofs == S->Mesh()->space_dimension()) {
    // diagonal tensor
    for (int lcv_dof=0; lcv_dof!=ndofs; ++lcv_dof) {
      for (int c=0; c!=ncells; ++c) {
        (*K_)[c](lcv_dof, lcv_dof) = (*perm)("cell",lcv_dof,c);
      }
    }
  } else {
    // ERROR -- unknown perm type
    ASSERT(0);
  }
};


// -----------------------------------------------------------------------------
// Update elemental discretization matrices with gravity terms.
//
// Must be called before applying boundary conditions and global assembling.
// -----------------------------------------------------------------------------
void Richards::AddGravityFluxes_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<Operators::MatrixMFD>& matrix) {

  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");

  // Get the rel perm, and ensure it is up to date.
  S->GetFieldEvaluator("relative_permeability")->HasFieldChanged(S.ptr(), "richards_pk");
  Teuchos::RCP<const CompositeVector> Krel = S->GetFieldData("numerical_rel_perm");

  // Get the density, in a mass basis, and ensure it is up to date.
  S->GetFieldEvaluator("mass_density_liquid")->HasFieldChanged(S.ptr(), "richards_pk");
  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");

  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int c_owned = rho->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    S->Mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
    double& Fc = matrix->Fc_cells()[c];

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->Mesh()->face_normal(f);

      double outward_flux = ( ((*K_)[c] * gravity) * normal) * dirs[n]
          * (*Krel)("face",f) *  (*Krel)("cell",c) * (*rho)("cell",c);
      Ff[n] += outward_flux;
      Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
    }
  }
};


// -----------------------------------------------------------------------------
// Updates global Darcy vector calculated by a discretization method.
// -----------------------------------------------------------------------------
void Richards::AddGravityFluxesToVector_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& darcy_flux) {

  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");

  // Get the rel perm, and ensure it is up to date.
  S->GetFieldEvaluator("relative_permeability")->HasFieldChanged(S.ptr(), "richards_pk");
  Teuchos::RCP<const CompositeVector> Krel = S->GetFieldData("numerical_rel_perm");

  // Get the density, in a mass basis, and ensure it is up to date.
  S->GetFieldEvaluator("mass_density_liquid")->HasFieldChanged(S.ptr(), "richards_pk");
  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");

  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int f_used = darcy_flux->size("face", true);
  int f_owned = darcy_flux->size("face", false);
  std::vector<bool> done(f_used, false);

  int c_owned = rho->size("cell");
  for (int c=0; c!=c_owned; ++c) {
    S->Mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->Mesh()->face_normal(f);

      if (f<f_owned && !done[f]) {
        (*darcy_flux)("face",f) += (((*K_)[c] * gravity) * normal)
            * (*Krel)("cell",c) * (*Krel)("face",f) * (*rho)("cell",c);
        done[f] = true;
      }
    }
  }
};

} //namespace
} //namespace
