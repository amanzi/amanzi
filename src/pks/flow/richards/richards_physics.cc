/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
  A base two-phase, thermal Richard's equation with water vapor.

  License: BSD
  Authors: Neil Carlson (version 1)
  Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
  Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */

#include "FieldEvaluator.hh"
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
  matrix_->CreateMFDrhsVectors();

  // derive fluxes
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");


  if (update_flux_ == UPDATE_FLUX_ITERATION) {
    // update the flux
    Teuchos::RCP<CompositeVector> flux =
        S->GetFieldData("darcy_flux", name_);

    matrix_->DeriveFlux(*pres, flux.ptr());
    if (matrix_->method() != Amanzi::Operators::FV_TPFA)
      AddGravityFluxesToVector_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), flux.ptr());
  }

  // assemble the stiffness matrix
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());

  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  // calculate the residual
  matrix_->ComputeNegativeResidual(*pres, g.ptr());

  Epetra_MultiVector sol_c = *(*pres).ViewComponent("cell");
    

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

  db_->WriteVector("res (acc)", g, true);
};


// ---------------------------------------------------------------------
// Add in mass source, in units of mol / m^3 s
// ---------------------------------------------------------------------
void Richards::AddSources_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  // external sources of energy
  if (is_source_term_) {
    Epetra_MultiVector& g_c = *g->ViewComponent("cell",false);

    // Update the source term
    S->GetFieldEvaluator("mass_source")->HasFieldChanged(S, name_);
    const Epetra_MultiVector& source1 =
        *S->GetFieldData("mass_source")->ViewComponent("cell",false);

    const Epetra_MultiVector& cv =
        *S->GetFieldData("cell_volume")->ViewComponent("cell",false);

    // Add into residual
    unsigned int ncells = g_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      g_c[0][c] -= source1[0][c] * cv[0][c];
    }

    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "Adding external source term" << std::endl;
      db_->WriteVector("  Q_ext", S->GetFieldData("mass_source").ptr(), false);
    }  
    db_->WriteVector("res (src)", g, false);
  }
}


void Richards::AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h) {
  // external sources of energy (temperature dependent source)
  if (is_source_term_ && !explicit_source_ &&
      S->GetFieldEvaluator("mass_source")->IsDependency(S, key_)) {
    std::vector<double>& Acc_cells = mfd_preconditioner_->Acc_cells();

    S->GetFieldEvaluator("mass_source")->HasFieldDerivativeChanged(S, name_, key_);
    const Epetra_MultiVector& dsource_dp =
        *S->GetFieldData("dmass_source_dpressure")->ViewComponent("cell",false);
    unsigned int ncells = dsource_dp.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      Acc_cells[c] -= dsource_dp[0][c];
    }
  }
}


// -------------------------------------------------------------
// Convert abs perm vector to tensor.
// -------------------------------------------------------------
void Richards::SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S) {
  // currently assumes isotropic perm, should be updated
  S->GetFieldEvaluator("permeability")->HasFieldChanged(S.ptr(), name_);
  const Epetra_MultiVector& perm = *S->GetFieldData("permeability")
      ->ViewComponent("cell",false);
  unsigned int ncells = perm.MyLength();
  unsigned int ndofs = perm.NumVectors();

  if (ndofs == 1) { // isotropic
    for (unsigned int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = perm[0][c] * perm_scale_;
    }
  } else if (ndofs == 2 && S->GetMesh()->space_dimension() == 3) {
    // horizontal and vertical perms
    for (unsigned int c=0; c!=ncells; ++c) {
      (*K_)[c](0, 0) = perm[0][c] * perm_scale_;
      (*K_)[c](1, 1) = perm[0][c] * perm_scale_;
      (*K_)[c](2, 2) = perm[1][c] * perm_scale_;
    }
  } else if (ndofs == S->GetMesh()->space_dimension()) {
    // diagonal tensor
    for (unsigned int lcv_dof=0; lcv_dof!=ndofs; ++lcv_dof) {
      for (unsigned int c=0; c!=ncells; ++c) {
        (*K_)[c](lcv_dof, lcv_dof) = perm[lcv_dof][c] * perm_scale_;
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

  if (matrix->method() == Amanzi::Operators::FV_TPFA){
    Teuchos::Ptr<Operators::Matrix_TPFA> matrix_tpfa = Teuchos::ptr_dynamic_cast<Operators::Matrix_TPFA>(matrix);
    AddGravityFluxes_FV_(g_vec, rel_perm, rho, matrix_tpfa);
    return;
  }




  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  if (rel_perm == Teuchos::null) { // no rel perm
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * rho_v[0][c];
        Ff[n] += outward_flux;
        Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
      }
    }

  } else if (!rel_perm->HasComponent("face")) { // rel perm on cells only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * krel_cells[0][c] * rho_v[0][c];
        Ff[n] += outward_flux;
        Fc -= outward_flux;  // Nonzero-sum contribution when not upwinding
      }
    }

  } else if (!rel_perm->HasComponent("cell")) { // rel perm on faces only
    rel_perm->ScatterMasterToGhosted("face");

    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * rho_v[0][c];
        Ff[n] += outward_flux * (scaled_constraint_ ? 1. : krel_faces[0][f]);
        Fc -= outward_flux * krel_faces[0][f] ;  // Nonzero-sum contribution when not upwinding
      }
    }

  } else { // rel perm on both cells and faces
    rel_perm->ScatterMasterToGhosted("face");

    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      double& Fc = matrix->Fc_cells()[c];

      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        double outward_flux = ( ((*K_)[c] * gravity) * normal) * krel_cells[0][c] * rho_v[0][c];

        //std::cout<<"grav "<<( ((*K_)[c] * gravity) * normal) * dirs[n]<<" rho "<<rho_v[0][c]<<" krel "<<krel_cells[0][c]<<"\n";

        Ff[n] += outward_flux * (scaled_constraint_ ? 1. : krel_faces[0][f]);
        Fc -= outward_flux * krel_faces[0][f];  // Nonzero-sum contribution when not upwinding

      }      
    }
  }
};

// -----------------------------------------------------------------------------
// Update elemental discretization matrices with gravity terms.
//
// Must be called before applying boundary conditions and global assembling.
// -----------------------------------------------------------------------------
void Richards::AddGravityFluxes_FV_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
        const Teuchos::Ptr<const CompositeVector>& rel_perm,
        const Teuchos::Ptr<const CompositeVector>& rho,
        const Teuchos::Ptr<Operators::Matrix_TPFA>& matrix) {

  int dim = g_vec->MyLength();
  AmanziGeometry::Point gravity(dim);
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  //Teuchos::RCP<const CompositeVector> rhs = matrix->rhs();
  //Epetra_MultiVector& F_cell = *rhs.ViewComponent("cell",false);

  std::vector<double>& Fc_cell = matrix->Fc_cells();

  Teuchos::RCP<Epetra_Vector> grav_terms = matrix->gravity_terms();

  AmanziMesh::Entity_ID_List faces, cells;
  std::vector<int> dirs;

  if (rel_perm == Teuchos::null) { // no rel perm

    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell", true);
    unsigned int ncells = rho->size("cell",false);

    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();
      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      Fc_cell[c] = 0.;

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
        double rho_avr = 0.;
        for (int i=0; i<cells.size(); ++i) rho_avr += rho_v[0][cells[i]];
        rho_avr *= 1./cells.size();

        double grav_flux = dirs[n] * gravity[dim-1] * (*grav_terms)[f] * rho_avr;
        if (cells.size() == 1){
          if (bc_markers_[f] != Amanzi::Operators::MATRIX_BC_DIRICHLET){
            Ff[n] -= grav_flux;
          }
          Fc_cell[c] -= grav_flux;
        }
        else{
          Fc_cell[c] -= grav_flux;
        }  
      }
    }

  }
  // else if (!rel_perm->HasComponent("cell")) { // rel perm on faces only
  else{
    rel_perm->ScatterMasterToGhosted("face");

    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell", true);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    unsigned int ncells = rho->size("cell",false);
    //Epetra_MultiVector& rhs_cells = *rhs_->ViewComponent("cell",false);

    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();
      Epetra_SerialDenseVector& Ff = matrix->Ff_cells()[c];
      Fc_cell[c] = 0.;


      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        mesh_->face_get_cells(f, AmanziMesh::USED, &cells);


        double rho_avr = 0.;
        for (int i=0; i<cells.size(); ++i) rho_avr += rho_v[0][cells[i]];
        rho_avr *= 1./cells.size();

        double grav_flux = dirs[n] * gravity[dim-1] * (*grav_terms)[f] * rho_avr * krel_faces[0][f];

        if (cells.size() == 1){
          if (bc_markers_[f] != Amanzi::Operators::MATRIX_BC_DIRICHLET){
            Ff[n] -= grav_flux;
          }
          Fc_cell[c] -= grav_flux;
        }
        else{
          Fc_cell[c] -= grav_flux;
        }  

      }
    }
  }
 // else {
 //   Errors::Message message(std::string("FV discretization doesn't support this type of relative permeability\n"));
 // }

}


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
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal) * dirs[n] * rho_v[0][c];
          done[f] = true;
        }
      }
    }

  } else if (!rel_perm->HasComponent("face")) { // rel perm on cells only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal) * dirs[n]
              * krel_cells[0][c] * rho_v[0][c];
          done[f] = true;
        }
      }
    }

  } else if (!rel_perm->HasComponent("cell")) { // rel perm on faces only
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal) * dirs[n]
              * krel_faces[0][f] * rho_v[0][c];
          done[f] = true;
        }
      }
    }

  } else { // rel perm on both cells and faces
    const Epetra_MultiVector& rho_v = *rho->ViewComponent("cell",false);
    const Epetra_MultiVector& krel_faces = *rel_perm->ViewComponent("face",true);
    const Epetra_MultiVector& krel_cells = *rel_perm->ViewComponent("cell",false);
    unsigned int ncells = rho->size("cell",false);
    for (unsigned int c=0; c!=ncells; ++c) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (unsigned int n=0; n!=faces.size(); ++n) {
        int f = faces[n];
        AmanziGeometry::Point normal = mesh_->face_normal(f) * dirs[n];
        if (tpfa_) {
          // normal must be vector connecting centroids, not true normal
          AmanziMesh::Entity_ID_List cells;
          mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
          const AmanziGeometry::Point& cell_centroid = mesh_->cell_centroid(c);
          AmanziGeometry::Point dc(mesh_->space_dimension());
          if (cells.size() == 1) {
            dc = mesh_->face_centroid(f) - cell_centroid;
          } else if (cells[0] == c) {
            dc = mesh_->cell_centroid(cells[1]) - cell_centroid;
          } else {
            dc = mesh_->cell_centroid(cells[0]) - cell_centroid;
          }
          normal = AmanziGeometry::norm(normal) / AmanziGeometry::norm(dc) * dc;
        }

        if (f<f_owned && !done[f]) {
          darcy_flux_v[0][f] += (((*K_)[c] * gravity) * normal) * dirs[n]
              * krel_cells[0][c] * krel_faces[0][f] * rho_v[0][c];
          done[f] = true;
        }
      }
    }
  }
};

// -------------------------------------------------------------
// Diffusion term, div -\phi s \tau n D grad \omega
// -------------------------------------------------------------
void Richards::AddVaporDiffusionResidual_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& g) {

  //res_vapor = Teuchos::rcp(new CompositeVector(*S->GetFieldData("pressure"))); 
  //res_vapor = Teuchos::rcp(new CompositeVector(*g)); 
  res_vapor->PutScalar(0.0);
  //Teuchos::RCP<CompositeVector> res_en = Teuchos::rcp(new CompositeVector(*g)); 
  //res_en->PutScalar(0.);

  // derive fluxes
  Teuchos::RCP<const CompositeVector> pres   = S->GetFieldData("pressure");
  Teuchos::RCP<const CompositeVector> temp   = S->GetFieldData("temperature");


  Teuchos::RCP<CompositeVector> vapor_diff_pres = S->GetFieldData("vapor_diffusion_pressure", name_);
  Teuchos::RCP<CompositeVector> vapor_diff_temp = S->GetFieldData("vapor_diffusion_temperature", name_);

  ///****** Compute contribution for pressure gradient

  //Epetra_MultiVector& coef_pr = *vapor_diff_pres->ViewComponent("cell",false);
  ComputeVaporDiffusionCoef(S, vapor_diff_pres, "pressure");

  // update the stiffness matrix
  matrix_vapor_->CreateMFDstiffnessMatrices(vapor_diff_pres.ptr());
  matrix_vapor_->CreateMFDrhsVectors();
  // assemble the stiffness matrix
  //matrix_vapor_->ApplyBoundaryConditions(bc_markers_, bc_values_, false);
  //  matrix_vapor_->AssembleGlobalMatrices();
  // calculate the residual
  matrix_vapor_->ComputeNegativeResidual(*pres, res_vapor.ptr());

  g->Update(1., *res_vapor, 1.);

  res_vapor->PutScalar(0.0);

  ///****** Compute contribution for temperature gradient
  Epetra_MultiVector& coef_tm = *S->GetFieldData("vapor_diffusion_temperature", name_)
                                   ->ViewComponent("cell",false);

  ComputeVaporDiffusionCoef(S, vapor_diff_temp, "temperature");
  // update the stiffness matrix
  matrix_vapor_->CreateMFDstiffnessMatrices(vapor_diff_temp.ptr());
  matrix_vapor_->CreateMFDrhsVectors();
  // assemble the stiffness matrix
  //matrix_vapor_->ApplyBoundaryConditions(bc_markers_, bc_values_, false);
  //  matrix_vapor_->AssembleGlobalMatrices();
  // calculate the residual
  matrix_vapor_->ComputeNegativeResidual(*temp, res_vapor.ptr());

  g->Update(1., *res_vapor, 1.);


}

  void Richards::ComputeVaporDiffusionCoef(const Teuchos::Ptr<State>& S, 
                                          Teuchos::RCP<CompositeVector>& vapor_diff, 
                                          std::string var_name){

   Epetra_MultiVector& diff_coef = *vapor_diff->ViewComponent("cell",false);

   S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S.ptr(), name_);
   const Epetra_MultiVector& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell",false);

   S->GetFieldEvaluator("molar_density_gas")->HasFieldChanged(S.ptr(), name_);
   const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell",false);

   S->GetFieldEvaluator("porosity")->HasFieldChanged(S.ptr(), name_);
   const Epetra_MultiVector& phi = *S->GetFieldData("porosity")->ViewComponent("cell",false);

   S->GetFieldEvaluator("saturation_gas")->HasFieldChanged(S.ptr(), name_);
   const Epetra_MultiVector& s_g = *S->GetFieldData("saturation_gas")->ViewComponent("cell",false);

   S->GetFieldEvaluator("mol_frac_gas")->HasFieldChanged(S.ptr(), name_);
   const Epetra_MultiVector& mlf_g = *S->GetFieldData("mol_frac_gas")->ViewComponent("cell",false);

   std::string key_t = "temperature";
   S->GetFieldEvaluator("mol_frac_gas")->HasFieldDerivativeChanged(S.ptr(), name_, key_t);
   const Epetra_MultiVector& dmlf_g_dt = *S->GetFieldData("dmol_frac_gas_dtemperature")->ViewComponent("cell",false);

   const Epetra_MultiVector& temp = *S->GetFieldData("temperature")->ViewComponent("cell",false);
   const Epetra_MultiVector& pressure = *S->GetFieldData("pressure")->ViewComponent("cell",false);
   const double& Patm = *S->GetScalarData("atmospheric_pressure");
   const double R = 8.3144621;

   unsigned int ncells = diff_coef.MyLength();

   const double a = 4./3.;
   const double b = 10./3.;
   const double D_ref = 0.282;
   const double P_ref = Patm;
   const double T_ref = 298;
   double D;

   for (unsigned int c=0; c!=ncells; ++c){

     D = D_ref*(P_ref/Patm)*pow(temp[0][c]/T_ref, 1.8);

     diff_coef[0][c] = D*pow(phi[0][c], a)*pow(s_g[0][c], b)*n_g[0][c];
     diff_coef[0][c] *= exp(-(Patm - pressure[0][c])/(n_l[0][c]*R*temp[0][c]));
   }

   if (var_name == "pressure"){
     //cout<<"Pressure vapor_diff\n";
     for (unsigned int c=0; c!=ncells; ++c){
       diff_coef[0][c] *= mlf_g[0][c] * (1./ (n_l[0][c]*R*temp[0][c]));
       //diff_coef[0][c] *= 0.;
       //cout<<diff_coef[0][c]<<" ";
     }
     //cout<<endl;
   }
   else if (var_name == "temperature"){
     for (unsigned int c=0; c!=ncells; ++c){
       diff_coef[0][c] *= (1./Patm)*dmlf_g_dt[0][c] + mlf_g[0][c]* (Patm - pressure[0][c])/ (n_l[0][c]*R*temp[0][c]*temp[0][c]);
       //diff_coef[0][c] =0;
     }
     
   }
   else{
     // Unknown variable name
     ASSERT(0);
   }    

}

} //namespace
} //namespace
