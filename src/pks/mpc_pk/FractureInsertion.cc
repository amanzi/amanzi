/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Insertion of fracture into a model changes mesh topology, mesh maps, etc.
  This class provides supporting tools and data structures that can be
  shared between MPC PKs.
*/

#include "CompositeVector.hh"
#include "UniqueLocalIndex.hh"

#include "FractureInsertion.hh"

namespace Amanzi {

/* *******************************************************************
* Constructor
******************************************************************* */
FractureInsertion::FractureInsertion(Teuchos::RCP<const AmanziMesh::Mesh>& mesh_matrix,
                                     Teuchos::RCP<const AmanziMesh::Mesh>& mesh_fracture)
  : mesh_matrix_(mesh_matrix), mesh_fracture_(mesh_fracture)
{}


/* *******************************************************************
* Inialization for matrix faces coupled to fracture cells.
******************************************************************* */
void
FractureInsertion::InitMatrixFaceToFractureCell(Teuchos::RCP<const Epetra_BlockMap> mmap,
                                                Teuchos::RCP<const Epetra_BlockMap> gmap)
{
  mmap_ = mmap;
  int npoints_owned = mmap_->NumMyPoints();

  cvs_matrix_ = Teuchos::rcp(new CompositeVectorSpace());
  cvs_fracture_ = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix_->SetMesh(mesh_matrix_)
    ->SetGhosted(true)
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, mmap, gmap, 1);

  cvs_fracture_->SetMesh(mesh_fracture_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  int ncells_owned_f =
    mesh_fracture_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  inds_matrix_ = std::make_shared<std::vector<std::vector<int>>>(npoints_owned);
  inds_fracture_ = std::make_shared<std::vector<std::vector<int>>>(npoints_owned);
  values_ = std::make_shared<std::vector<double>>(npoints_owned);

  int np(0);
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    int first = mmap_->FirstPointInElement(f);
    int ndofs = mmap_->ElementSize(f);

    for (int k = 0; k < ndofs; ++k) {
      (*inds_matrix_)[np].resize(1);
      (*inds_fracture_)[np].resize(1);
      (*inds_matrix_)[np][0] = first + k;
      (*inds_fracture_)[np][0] = c;
      np++;
    }
  }

  inds_matrix_->resize(np);
  inds_fracture_->resize(np);
  values_->resize(np);
}


/* *******************************************************************
* Inialization matrix cell coupled to fracture cells.
******************************************************************* */
void
FractureInsertion::InitMatrixCellToFractureCell()
{
  cvs_matrix_ = Teuchos::rcp(new CompositeVectorSpace());
  cvs_fracture_ = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix_->SetMesh(mesh_matrix_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  cvs_fracture_->SetMesh(mesh_fracture_)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // -- indices are fluxes on matrix-fracture interface
  int ncells_owned_f =
    mesh_fracture_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  inds_matrix_ = std::make_shared<std::vector<std::vector<int>>>(2 * ncells_owned_f);
  inds_fracture_ = std::make_shared<std::vector<std::vector<int>>>(2 * ncells_owned_f);
  values_ = std::make_shared<std::vector<double>>(2 * ncells_owned_f, 0.0);

  int np(0), ierr(0);

  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    auto cells = mesh_matrix_->getFaceCells(f);
    int ncells = cells.size();
    if (ncells != 2) ierr = 1;

    for (int k = 0; k < ncells; ++k) {
      (*inds_matrix_)[np].resize(1);
      (*inds_fracture_)[np].resize(1);
      (*inds_matrix_)[np][0] = cells[k];
      (*inds_fracture_)[np][0] = c;
      np++;
    }
  }

  // if at least one process errors, all must terminate
  int recv = 0;
  mesh_matrix_->getComm()->MaxAll(&ierr, &recv, 1);
  if (recv != 0) {
    Errors::Message msg("A fracture face has only one attached cell, insteaf of two.");
    Exceptions::amanzi_throw(msg);
  }
}


/* *******************************************************************
* Set up data
******************************************************************* */
void
FractureInsertion::SetValues(const Epetra_MultiVector& kn, double scale)
{
  int np(0);
  int ncells_owned_f =
    mesh_fracture_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  int kmax = kn.NumVectors();

  for (int c = 0; c < ncells_owned_f; ++c) {
    double area = mesh_fracture_->getCellVolume(c);
    int f = mesh_fracture_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    int ndofs = mmap_->ElementSize(f);

    for (int k = 0; k < ndofs; ++k) {
      int k1 = std::min(k, kmax - 1);
      (*values_)[np] = kn[k1][c] * area * scale;
      np++;
    }
  }
}


/* *******************************************************************
* Coupling coefficients depend on flux sign
******************************************************************* */
void
FractureInsertion::SetValues(const CompositeVector& flux)
{
  int np(0), dir;
  int ncells_owned_f =
    mesh_fracture_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  if (!values2_.get()) {
    values2_ = std::make_shared<std::vector<double>>(2 * ncells_owned_f, 0.0);
  }

  const auto& flux_f = *flux.ViewComponent("face");
  const auto& mmap = flux.Map().Map("face", false);

  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    int first = mmap->FirstPointInElement(f);

    auto cells = mesh_matrix_->getFaceCells(f);
    int ncells = cells.size();
    mesh_matrix_->getFaceNormal(f, cells[0], &dir);
    int shift = Operators::UniqueIndexFaceToCells(*mesh_matrix_, f, cells[0]);

    for (int k = 0; k < ncells; ++k) {
      // since cells could be ordered differenty then points, we need a map
      double tmp = flux_f[0][first + shift] * dir;

      if (tmp > 0) {
        (*values_)[np] = tmp;
        (*values2_)[np] = 0.0;
      } else {
        (*values_)[np] = 0.0;
        (*values2_)[np] = -tmp;
      }

      dir = -dir;
      shift = 1 - shift;
      np++;
    }
  }
}

} // namespace Amanzi
