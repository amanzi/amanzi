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
  int ncells_f = fractures->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_m = mstk->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  std::vector<int> points(nfaces_m, 1);

  for (int c = 0; c < ncells_f; ++c) {
    int f = fractures->entity_get_parent(AmanziMesh::CELL, c);
    mstk->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    // points[f] = cells.size();
  }

  // create master map with two points on each fracture face
  auto& mfmap = mesh_->face_map(false);
  int nglobal = mfmap.NumGlobalElements();
  int nlocal = mfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  mfmap.MyGlobalElements(&gids[0]);

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(nglobal, nlocal, &gids[0], &points[0], 0, mfmap.Comm()));

  // create ghosted map with two points on each fracture face
  auto& gfmap = mesh_->face_map(true);
  nglobal = gfmap.NumGlobalElements();
  nlocal = gfmap.NumMyElements();

  points.resize(nlocal);
  gids.resize(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  for (int i = 0; i < nlocal; ++i) points[i] = 1;
  auto gmap = Teuchos::rcp(new Epetra_BlockMap(nglobal, nlocal, &gids[0], &points[0], 0, gfmap.Comm()));

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

}  // namespace Operators
}  // namespace Amanzi

