/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  Insertion of fracture into a model changes mesh topology, mesh maps, etc.
  This class provides supporting tools and data structures that can be
  shared between MPC PKs.
*/

#include "FractureInsertion.hh"

namespace Amanzi {

/* *******************************************************************
* Constructor
******************************************************************* */
FractureInsertion::FractureInsertion(
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh_matrix,
    Teuchos::RCP<const AmanziMesh::Mesh>& mesh_fracture)
  :  mesh_matrix_(mesh_matrix),
     mesh_fracture_(mesh_fracture)
{
}


/* *******************************************************************
* Inialization with optional block maps for new matrix mesh topology.
******************************************************************* */
void FractureInsertion::Init(
   Teuchos::RCP<const Epetra_BlockMap> mmap,
   Teuchos::RCP<const Epetra_BlockMap> gmap)
{
  mmap_ = mmap;
  int npoints_owned = mmap_->NumMyPoints();

  cvs_matrix_ = Teuchos::rcp(new CompositeVectorSpace());
  cvs_fracture_ = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix_->SetMesh(mesh_matrix_)->SetGhosted(true)
             ->AddComponent("face", AmanziMesh::FACE, mmap, gmap, 1);

  cvs_fracture_->SetMesh(mesh_fracture_)->SetGhosted(true)
               ->AddComponent("cell", AmanziMesh::CELL, 1);

  int ncells_owned_f = mesh_fracture_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  inds_matrix_ = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  inds_fracture_ = std::make_shared<std::vector<std::vector<int> > >(npoints_owned);
  values_ = std::make_shared<std::vector<double> >(npoints_owned);

  int np(0);
  for (int c = 0; c < ncells_owned_f; ++c) {
    int f = mesh_fracture_->entity_get_parent(AmanziMesh::CELL, c);
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
* Set up data
******************************************************************* */
void FractureInsertion::SetValues(const Epetra_MultiVector& kn, double scale)
{
  int np(0);
  int ncells_owned_f = mesh_fracture_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned_f; ++c) {
    double area = mesh_fracture_->cell_volume(c);
    int f = mesh_fracture_->entity_get_parent(AmanziMesh::CELL, c);
    int ndofs = mmap_->ElementSize(f);

    for (int k = 0; k < ndofs; ++k) {
      (*values_)[np] = kn[0][c] * area * scale;
      np++;
    }
  }
}

}  // namespace Amanzi
