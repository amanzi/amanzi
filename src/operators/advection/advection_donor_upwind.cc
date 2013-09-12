/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Donor upwind advection.
   ------------------------------------------------------------------------- */

#include "advection_donor_upwind.hh"

namespace Amanzi {
namespace Operators {


AdvectionDonorUpwind::AdvectionDonorUpwind(Teuchos::ParameterList& advect_plist,
        const Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    Advection(advect_plist, mesh) {

  f_begin_ = mesh_->face_map(true).MinLID();
  f_count_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  f_owned_ = f_begin_ + f_count_;
  f_end_ = mesh_->face_map(true).MaxLID() + 1;

  c_begin_ = mesh_->cell_map(true).MinLID();
  c_count_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  c_owned_ = c_begin_ + c_count_;
  c_end_ = mesh_->cell_map(true).MaxLID() + 1;

  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));

};


// set flux and determine upwind cells
void AdvectionDonorUpwind::set_flux(const Teuchos::RCP<const CompositeVector>& flux) {
  flux_ = flux;
  IdentifyUpwindCells_();
};


void AdvectionDonorUpwind::Apply(const Teuchos::RCP<Functions::BoundaryFunction>& bc_flux,
                                 bool include_bc_fluxes) {

  field_->ScatterMasterToGhosted("cell"); // communicate the cells

  // collect fluxes in faces
  Epetra_MultiVector& field_f = *field_->ViewComponent("face",true);
  Epetra_MultiVector& field_c = *field_->ViewComponent("cell", true);

  Teuchos::RCP<const CompositeVector> flux_const(flux_);
  const Epetra_MultiVector& flux = *flux_const->ViewComponent("face",false);

  for (int f=f_begin_; f != f_end_; ++f) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];

    if (c1 >=0) {
      double u = std::abs(flux[0][f]);
      for (int i=0; i != num_dofs_; ++i) {
        field_f[i][f] = u * field_c[i][c1];
      }
    }
  }

  // patch up Neumann bcs -- only works for 1 dof?
  for (Functions::BoundaryFunction::Iterator bc = bc_flux->begin();
       bc!=bc_flux->end(); ++bc) {
    if (include_bc_fluxes) {
      if ((*upwind_cell_)[bc->first] >= 0) {
        field_f[0][bc->first] = bc->second*mesh_->face_area(bc->second);
      } else {
        field_f[0][bc->first] = -bc->second*mesh_->face_area(bc->second);
      }
    } else {
      // HACKED for fluxes included in diffusion term.  This needs
      // rethought. --etc
      field_f[0][bc->first] = 0.;
    }
  }

  field_c.PutScalar(0.);
  // put fluxes in cell
  for (int f=f_begin_; f!=f_end_; ++f) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    if (c1 >=0 && c1 < c_owned_) {
      for (int i=0; i<num_dofs_; i++) {
        field_c[i][c1] -= field_f[i][f];
      }
    }
    if (c2 >=0 && c2 < c_owned_) {
      for (int i=0; i<num_dofs_; i++) {
        field_c[i][c2] += field_f[i][f];
      }
    }
  }
};


void AdvectionDonorUpwind::IdentifyUpwindCells_() {

  upwind_cell_->PutValue(-1);
  downwind_cell_->PutValue(-1);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;
  flux_->ScatterMasterToGhosted("face");
  Teuchos::RCP<const CompositeVector> flux_const(flux_);
  const Epetra_MultiVector& flux_f = *flux_const->ViewComponent("face",true);

  for (int c=c_begin_; c != c_end_; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (unsigned int i = 0; i != faces.size(); ++i) {
      AmanziMesh::Entity_ID f = faces[i];
      if (flux_f[0][f] * fdirs[i] >= 0) {
        (*upwind_cell_)[f] = c;
      } else {
        (*downwind_cell_)[f] = c;
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi
