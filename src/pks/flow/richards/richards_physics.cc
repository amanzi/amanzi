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
// Diffusion term, div K grad (p + rho*g*z)
// -------------------------------------------------------------
void Richards::ApplyDiffusion_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  // update the rel perm according to the scheme of choice
  bool update = UpdatePermeabilityData_(S.ptr());

  // update the stiffness matrix
  Teuchos::RCP<const CompositeVector> rel_perm =
    S->GetFieldData("numerical_rel_perm");
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());

  // derive fluxes
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
  if (update_flux_ == UPDATE_FLUX_ITERATION ||
      coupled_to_surface_via_head_ || coupled_to_surface_via_full_) {
    // update the flux
    Teuchos::RCP<CompositeVector> flux =
        S->GetFieldData("darcy_flux", name_);

    // derive flux
    if (coupled_to_surface_via_head_) {
      // this deals with the lag associated with teh coupling?
      Teuchos::RCP<CompositeVector> pres_corrected =
          Teuchos::rcp(new CompositeVector(*pres));
      *pres_corrected = *pres;

      Epetra_MultiVector& pres_corrected_f = *pres_corrected->ViewComponent("face",false);
      const Epetra_MultiVector& head = *S->GetFieldData("surface_pressure")
          ->ViewComponent("cell",false);
      Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh("surface");

      int ncells_surface = head.MyLength();
      for (int c=0; c!=ncells_surface; ++c) {
        // -- get the surface cell's equivalent subsurface face and neighboring cell
        AmanziMesh::Entity_ID f =
            surface->entity_get_parent(AmanziMesh::CELL, c);
        pres_corrected_f[0][f] = head[0][c];
      }

      matrix_->DeriveFlux(*pres_corrected, flux.ptr());
    } else {
      matrix_->DeriveFlux(*pres, flux.ptr());
    }

    AddGravityFluxesToVector_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), flux.ptr());
    flux->ScatterMasterToGhosted();
  }

  // assemble the stiffness matrix
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres, g.ptr());
};


// -------------------------------------------------------------
// Accumulation of water term du/dt
// -------------------------------------------------------------
void Richards::AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) {
  double dt = S_next_->time() - S_inter_->time();

  // update the water content at both the old and new times.
  S_next_->GetFieldEvaluator("water_content")->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator("water_content")->HasFieldChanged(S_inter_.ptr(), name_);

  // get these fields
  Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
  Teuchos::RCP<const CompositeVector> wc0 = S_inter_->GetFieldData("water_content");

  // Water content only has cells, while the residual has cells and faces.
  g->ViewComponent("cell",false)->Update(1.0/dt, *wc1->ViewComponent("cell",false),
          -1.0/dt, *wc0->ViewComponent("cell",false), 1.0);
};


// -------------------------------------------------------------
// Convert abs perm vector to tensor.
// -------------------------------------------------------------
void Richards::SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S) {
  // currently assumes isotropic perm, should be updated
  S->GetFieldEvaluator("permeability")->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> perm = S->GetFieldData("permeability");
  int ncells = perm->size("cell");
  int ndofs = perm->num_dofs("cell");

  if (ndofs == 1) { // isotropic
    for (int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = (*perm)("cell",c);
    }
  } else if (ndofs == 2 && S->GetMesh()->space_dimension() == 3) {
    // horizontal and vertical perms
    for (int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = (*perm)("cell",0,c);
      (*K_)[c](1, 1) = (*perm)("cell",0,c);
      (*K_)[c](2, 2) = (*perm)("cell",1,c);
    }
  } else if (ndofs == S->GetMesh()->space_dimension()) {
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
void Richards::AddGravityFluxes_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
        const Teuchos::Ptr<const CompositeVector>& rel_perm,
        const Teuchos::Ptr<const CompositeVector>& rho,
        const Teuchos::Ptr<Operators::MatrixMFD>& matrix) {

  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  if (rel_perm == Teuchos::null) { // no rel perm
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * dirs[n]
            * rho_v[0][c];
        Ff[n] += outward_flux;
        Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
      }
    }

  } else if (!rel_perm->has_component("face")) { // rel perm on cells only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      rho->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = rho->mesh()->face_normal(f);

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * dirs[n]
            * krel_cells[0][c] * rho_v[0][c];
        Ff[n] += outward_flux;
        Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
      }
    }

  } else if (!rel_perm->has_component("cell")) { // rel perm on faces only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      rho->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = rho->mesh()->face_normal(f);

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * dirs[n]
            * krel_faces[0][f] * rho_v[0][c];
        Ff[n] += outward_flux;
        Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
      }
    }

  } else { // rel perm on both cells and faces
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      rho->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = rho->mesh()->face_normal(f);

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * dirs[n]
            * krel_faces[0][f] * krel_cells[0][c] * rho_v[0][c];
        Ff[n] += outward_flux;
        Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
      }
    }
  }
};


// -----------------------------------------------------------------------------
// Updates global Darcy vector calculated by a discretization method.
// -----------------------------------------------------------------------------
void Richards::AddGravityFluxesToVector_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
        const Teuchos::Ptr<const CompositeVector>& rel_perm,
        const Teuchos::Ptr<const CompositeVector>& rho,
        const Teuchos::Ptr<CompositeVector>& darcy_flux) {

  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int f_owned = darcy_flux->size("face", false);
  std::vector<bool> done(darcy_flux->size("face",true), false);
  Epetra_MultiVector& darcy_flux_v = *darcy_flux->ViewComponent("face",false);

  if (rel_perm == Teuchos::null) { // no rel perm
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal) * rho_v[0][c];
          done[f] = true;
        }
      }
    }

  } else if (!rel_perm->has_component("face")) { // rel perm on cells only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      darcy_flux->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = darcy_flux->mesh()->face_normal(f);
        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal)
              * krel_cells[0][c] * rho_v[0][c];
          done[f] = true;
        }
      }
    }

  } else if (!rel_perm->has_component("cell")) { // rel perm on faces only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      darcy_flux->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = darcy_flux->mesh()->face_normal(f);
        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal)
              * krel_faces[0][f] * rho_v[0][c];
          done[f] = true;
        }
      }
    }

  } else { // rel perm on both cells and faces
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    int ncells = rho->size("cell",false);
    for (int c=0; c!=ncells; ++c) {
      darcy_flux->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& normal = darcy_flux->mesh()->face_normal(f);
        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal)
              * krel_cells[0][c] * krel_faces[0][f] * rho_v[0][c];
          done[f] = true;
        }
      }
    }
  }
};

} //namespace
} //namespace
