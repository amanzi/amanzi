/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_BlockMap.h"

// Amanzi
#include "Mesh_MSTK.hh"
#include "Op_Cell_FaceCell.hh"
#include "Operator_FaceCell.hh"
#include "PDE_DiffusionFracturedMatrix.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void PDE_DiffusionFracturedMatrix::Init(Teuchos::ParameterList& plist)
{
  // extract mesh in fractures
  std::vector<std::string> names = plist.get<Teuchos::Array<std::string> >("fractures").toVector();

  Teuchos::RCP<const AmanziMesh::Mesh_MSTK> mstk =
      Teuchos::rcp_static_cast<const AmanziMesh::Mesh_MSTK>(mesh_);
  Teuchos::RCP<const AmanziMesh::Mesh> fractures =
      Teuchos::rcp(new AmanziMesh::Mesh_MSTK(*mstk, names, AmanziMesh::FACE));

  // create list of parents for owned cells
  AmanziMesh::Entity_ID_List cells;
  int ncells_f = fractures->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_m = mstk->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  points_.resize(nfaces_m, 1);

  for (int c = 0; c < ncells_f; ++c) {
    int f = fractures->entity_get_parent(AmanziMesh::CELL, c);
    mstk->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    points_[f] = cells.size();
  }

  // create ghosted map with two points on each fracture face
  auto& gfmap = mesh_->face_map(true);
  int nlocal = gfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  auto gmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], &points_[0], 0, gfmap.Comm()));

  // create master map with two points on each fracture face
  auto& mfmap = mesh_->face_map(false);
  nlocal = mfmap.NumMyElements();

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], &points_[0], 0, mfmap.Comm()));

  // create global operator
  std::string compname("face");
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true);
  cvs->AddComponent(compname, mmap, gmap, 1);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);

  global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs, plist));

  std::string name = "DiffusionFracturedMatrix: CELL_FACE+CELL";
  local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
  global_op_->OpPushBack(local_op_);

  // other parameters
  scaled_constraint_ = false;
  little_k_ = OPERATOR_LITTLE_K_NONE;
  newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  exclude_primary_terms_ = false;
  mass_matrices_initialized_ = false;
  K_ = Teuchos::null;
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void PDE_DiffusionFracturedMatrix::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& flux,
    const Teuchos::Ptr<const CompositeVector>& u)
{
  PDE_DiffusionMFD::UpdateMatrices(flux, u);

  AmanziMesh::Entity_ID_List cells, faces;
  const auto& cmap = mesh_->cell_map(true);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    int npoints(0);
    std::vector<int> map;

    for (int i = 0; i < nfaces; ++i) {
      int f = faces[i];
      int n = points_[f];

      int shift(0);
      if (n == 2) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        if (cells.size() == 2) {
          int gid = cmap.GID(c);
          int gid_min = std::min(cmap.GID(cells[0]), cmap.GID(cells[1]));
          if (gid > gid_min) shift = 1;
        }
      }

      map.push_back(npoints + shift);
      npoints += n;
    }
    
    // resize matrix
    if (npoints != nfaces) {
      auto& Acell = local_op_->matrices[c];

      WhetStone::DenseMatrix Anew(npoints, npoints);
      Anew.PutScalar(0.0);

      for (int i = 0; i < nfaces; ++i) {
        for (int j = 0; j < nfaces; ++j) {
          Anew(map[i], map[j]) = Acell(i, j);
        }
      }

      local_op_->matrices[c] = Anew;
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi

