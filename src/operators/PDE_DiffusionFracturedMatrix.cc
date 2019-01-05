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
#include "CompositeVector.hh"
#include "Mesh_MSTK.hh"
#include "Op_Cell_FaceCell.hh"
#include "Operator_FaceCell.hh"
#include "ParallelCommunication.hh"
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

  points_ = Teuchos::rcp(new Epetra_IntVector(mesh_->face_map(true)));
  points_->PutValue(1);

  for (int c = 0; c < ncells_f; ++c) {
    int f = fractures->entity_get_parent(AmanziMesh::CELL, c);
    mstk->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    (*points_)[f] = cells.size();
  }

  ParallelCommunication pp(mesh_);
  pp.CopyMasterFace2GhostFace(*points_);

  // create ghosted map with two points on each fracture face
  auto& gfmap = mesh_->face_map(true);
  int nlocal = gfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  int* data; 
  points_->ExtractView(&data);
  auto gmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, gfmap.Comm()));

  // create master map with two points on each fracture face
  auto& mfmap = mesh_->face_map(false);
  nlocal = mfmap.NumMyElements();

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, mfmap.Comm()));

  // create global operator
  std::string compname("face");
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true);
  cvs->AddComponent(compname, mmap, gmap, 1);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);

  global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs, plist));
  global_op_->set_variable_dofs(true);

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
      int n = (*points_)[f];

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
    map.push_back(npoints);
    
    // resize matrix
    if (npoints != nfaces) {
      auto& Acell = local_op_->matrices[c];

      WhetStone::DenseMatrix Anew(npoints + 1, npoints + 1);
      Anew.PutScalar(0.0);

      for (int i = 0; i < nfaces + 1; ++i) {
        for (int j = 0; j < nfaces + 1; ++j) {
          Anew(map[i], map[j]) = Acell(i, j);
        }
      }

      local_op_->matrices[c] = Anew;
    }
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: eliminate them or not. Finally we may add the essential BC
* the the system of equations as the trivial equations.
*
* Note: BCs imposed incorrectly on faces with many DOFs.
****************************************************************** */
void PDE_DiffusionFracturedMatrix::ApplyBCs(
    bool primary, bool eliminate, bool essential_eqn)
{
  // apply diffusion type BCs to FACE-CELL system
  AmanziMesh::Entity_ID_List faces;

  const std::vector<int>& bc_model_trial = bcs_trial_[0]->bc_model();
  const std::vector<int>& bc_model_test = bcs_test_[0]->bc_model();

  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  const std::vector<double>& bc_mixed = bcs_trial_[0]->bc_mixed();

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

  const auto& map = rhs_face.Map();

  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();
    
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    int nrows = Acell.NumRows();
        
    // essential conditions for test functions
    int np(0);
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        for (int k = 0; k < map.ElementSize(f); ++k) {
          for (int m = 0; m < nrows; m++) Acell(np + k, m) = 0.0;
        }
      }
      np += map.ElementSize(f);
    }

    // conditions for trial functions
    np = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int first = map.FirstPointInElement(f);
      double value = bc_value[f];

      if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) {
        // make a copy of elemental matrix for post-processing
        if (flag) {
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }

        if (eliminate) { 
          int mp(0);
          for (int m = 0; m < nfaces; m++) {
            int g = faces[m];
            int k = map.FirstPointInElement(g);
            rhs_face[0][k] -= Acell(mp, np) * value;
            Acell(mp, np) = 0.0;
            mp += map.ElementSize(g);
          }

          rhs_cell[0][c] -= Acell(nrows - 1, np) * value;
          Acell(nrows - 1, np) = 0.0;
        }

        if (essential_eqn) {
          rhs_face[0][first] = value;
          Acell(np, np) = 1.0;
        }

      } else if (bc_model_trial[f] == OPERATOR_BC_NEUMANN && primary) {
        rhs_face[0][first] -= value * mesh_->face_area(f);

      } else if (bc_model_trial[f] == OPERATOR_BC_TOTAL_FLUX && primary) {
        rhs_face[0][first] -= value * mesh_->face_area(f);

      } else if (bc_model_trial[f] == OPERATOR_BC_MIXED && primary) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        double area = mesh_->face_area(f);
        rhs_face[0][first] -= value * area;
        Acell(np, np) += bc_mixed[f] * area;
      }
      np += map.ElementSize(f);
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}

}  // namespace Operators
}  // namespace Amanzi

