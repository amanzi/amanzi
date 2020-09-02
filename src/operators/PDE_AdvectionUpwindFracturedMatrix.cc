/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#include <vector>

#include "WhetStoneDefs.hh"

#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "Op_Face_Cell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "PDE_AdvectionUpwindFracturedMatrix.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Advection requires a velocity field.
****************************************************************** */
void PDE_AdvectionUpwindFracturedMatrix::Setup(const CompositeVector& u)
{
  IdentifyUpwindCells_(u);
}

  
/* ******************************************************************
* A simple first-order transport method.
* Advection operator is of the form: div (u C), where u is the given
* velocity field and C is the advected field.
****************************************************************** */
void PDE_AdvectionUpwindFracturedMatrix::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& u)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  AmanziMesh::Entity_ID_List cells;
  const Epetra_MultiVector& uf = *u->ViewComponent("face");
  const auto& gmap = uf.Map();

  for (int f = 0; f < nfaces_owned; ++f) {
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    int g = gmap.FirstPointInElement(f);
    double umod = fabs(uf[0][g]);
    if (c1 < 0) {
      Aface(0, 0) = umod;
    } else if (c2 < 0) {
      Aface(0, 0) = umod;
    } else {
      int i = (cells[0] == c1) ? 0 : 1;
      Aface(i, i) = umod;
      Aface(1 - i, i) = -umod;
    }

    matrix[f] = Aface;
  }

  // removed matrices fof faces where fracture is located
  AmanziMesh::Entity_ID_List block;
  std::vector<double> vofs;
  for (int i = 0; i < fractures_.size(); ++i) {
    mesh_->get_set_entities_and_vofs(fractures_[i], AmanziMesh::FACE, 
                                     AmanziMesh::Parallel_type::OWNED, &block, &vofs);

    for (int n = 0; n < block.size(); ++n) {
      matrix[block[n]] *= 0.0;
    }
  }
}


/* ******************************************************************
* Add a simple first-order upwind method where the advected quantity
* is not the primary variable (used in Jacobians).
* Advection operator is of the form: div (q H(u))
*     q:    flux
*     H(u): advected quantity (i.e. enthalpy)
****************************************************************** */
void PDE_AdvectionUpwindFracturedMatrix::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& dhdT)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;

  AmanziMesh::Entity_ID_List cells;
  const Epetra_MultiVector& uf = *u->ViewComponent("face");
  const auto& gmap = uf.Map();

  dhdT->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& dh = *dhdT->ViewComponent("cell", true);

  for (int f = 0; f < nfaces_owned; ++f) {
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    WhetStone::DenseMatrix Aface(ncells, ncells);
    Aface.PutScalar(0.0);

    int g = gmap.FirstPointInElement(f);
    double umod = fabs(uf[0][g]);
    if (c1 < 0) {
      Aface(0, 0) = umod * dh[0][c2];
    } else if (c2 < 0) {
      Aface(0, 0) = umod * dh[0][c1];
    } else {
      int i = (cells[0] == c1) ? 0 : 1;
      Aface(i, i) = umod * dh[0][c1];
      Aface(1 - i, i) = -umod * dh[0][c2];
    }

    matrix[f] = Aface;
  }

  // removed matrices fof faces where fracture is located
  AmanziMesh::Entity_ID_List block;
  std::vector<double> vofs;
  for (int i = 0; i < fractures_.size(); ++i) {
    mesh_->get_set_entities_and_vofs(fractures_[i], AmanziMesh::FACE, 
                                     AmanziMesh::Parallel_type::OWNED, &block, &vofs);

    for (int n = 0; n < block.size(); ++n) {
      matrix[block[n]] *= 0.0;
    }
  }
}


/* *******************************************************************
* Identify flux direction based on orientation of the face normal 
* and sign of the  Darcy velocity.                               
******************************************************************* */
void PDE_AdvectionUpwindFracturedMatrix::IdentifyUpwindCells_(const CompositeVector& u)
{
  u.ScatterMasterToGhosted("face");
  const Epetra_MultiVector& uf = *u.ViewComponent("face", true);
  const auto& gmap = uf.Map();

  const Epetra_Map& fmap_wghost = mesh_->face_map(true);
  upwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell_ = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell_)[f] = -1;  // negative value indicates boundary
    (*downwind_cell_)[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      int g = gmap.FirstPointInElement(f);

      if (uf[0][g] * fdirs[i] >= 0) {
        (*upwind_cell_)[f] = c;
      } else {
        (*downwind_cell_)[f] = c;
      }
    }
  }
}


/* ******************************************************************
* Initialize additional parameters
****************************************************************** */
void PDE_AdvectionUpwindFracturedMatrix::InitAdvection_(Teuchos::ParameterList& plist)
{
  fractures_ = plist.get<Teuchos::Array<std::string> >("fracture").toVector();
}

}  // namespace Operators
}  // namespace Amanzi
