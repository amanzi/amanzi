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
        Teuchos::RCP<AmanziMesh::Mesh> mesh) :
    Advection(advect_plist, mesh) {

  f_begin_ = mesh_->face_map(true).MinLID();
  f_count_ = mesh_->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  f_owned_ = f_begin_ + f_count_;
  f_end_ = mesh_->face_map(true).MaxLID() + 1;

  c_begin_ = mesh_->cell_map(true).MinLID();
  c_count_ = mesh_->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  c_owned_ = c_begin_ + c_count_;
  c_end_ = mesh_->cell_map(true).MaxLID() + 1;

  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));

};


// set flux and determine upwind cells
void AdvectionDonorUpwind::set_flux(Teuchos::RCP<const CompositeVector>& flux) {
  flux_ = flux;
  IdentifyUpwindCells_();
};


void AdvectionDonorUpwind::Advect() {

  field_->ScatterMasterToGhosted("cell"); // communicate the cells

  double u;
  // collect fluxes in faces
  field_->ViewComponent("face")->PutScalar(0.0);
  for (int f=f_begin_; f != f_end_; ++f) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    if (c1 >=0) {
      u = fabs((*flux_)(f));
      for (int i=0; i != num_dofs_; ++i) {
        (*field_)("face",i,f) = u * (*field_)("cell",i,c1);
      }
    }
  }

  field_->ViewComponent("cell")->PutScalar(0.0);
  // put fluxes in cell
  for (int f=f_begin_; f != f_end_; ++f) {  // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    if (c1 >=0 && c1 < c_owned_ && c2 >= 0 && c2 < c_owned_) {
      for (int i=0; i<num_dofs_; i++) {
        double flux = (*field_)("face",i,f);
        (*field_)("cell",i,c1) -= flux;
        (*field_)("cell",i,c2) += flux;
      }
    } else if (c1 >=0 && c1 < c_owned_ && (c2 >= c_owned_ || c2 < 0)) {
      for (int i=0; i<num_dofs_; i++) {
        (*field_)("cell",i,c1) -= (*field_)("face",i,f);
      }
    } else if (c1 >= c_owned_ && c2 >= 0 && c2 < c_owned_) {
      for (int i=0; i<num_dofs_; i++) {
        (*field_)("cell",i,c2) += (*field_)("face",i,f);
      }
    }
  }
};


void AdvectionDonorUpwind::IdentifyUpwindCells_() {

  upwind_cell_->PutValue(-1);
  downwind_cell_->PutValue(-1);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c=c_begin_; c != c_end_; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i != faces.size(); ++i) {
      int f = faces[i];
      if ((*flux_)(f) * fdirs[i] >= 0) {
        (*upwind_cell_)[f] = c;
      } else {
        (*downwind_cell_)[f] = c;
      }
    }
  }
};

} // namespace Operators
} // namespace Amanzi
