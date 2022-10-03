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
#include "MeshFactory.hh"
#include "Op_Cell_FaceCell.hh"
#include "Operator_FaceCell.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void PDE_DiffusionFracturedMatrix::Init(Teuchos::ParameterList& plist)
{
  // extract mesh in fractures
  AmanziMesh::MeshFactory meshfactory(mesh_->get_comm(), mesh_->geometric_model());
  meshfactory.set_preference(AmanziMesh::Preference({AmanziMesh::Framework::MSTK}));

  std::vector<std::string> names = plist.get<Teuchos::Array<std::string> >("fracture").toVector();
  Teuchos::RCP<const AmanziMesh::Mesh> fracture = meshfactory.create(mesh_, names, AmanziMesh::FACE);

  // create global operator
  cvs_ = CreateFracturedMatrixCVS(mesh_, fracture);

  global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs_, plist));

  std::string name("DiffusionFracturedMatrix: CELL_FACE+CELL");
  local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
  global_op_->OpPushBack(local_op_);

  // little-k options
  name = plist.get<std::string>("nonlinear coefficient", "none");
  little_k_ = OPERATOR_LITTLE_K_NONE;
  if (name == "standard: cell") {
    little_k_ = OPERATOR_LITTLE_K_STANDARD;  // cell-centered scheme.
  } else if (name != "none") {
    Errors::Message msg;
    msg << "Diffusion solver does not support coefficient \"" << name << "\".\n";
    Exceptions::amanzi_throw(msg);
  }

  // other parameters
  gravity_ = plist.get<bool>("gravity");
  scaled_constraint_ = false;
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

  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  if (k_.get() && k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");

  const auto& fmap = *cvs_->Map("face", true);

  int dim = mesh_->space_dimension();
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);

  global_op_->rhs()->PutScalarGhosted(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    std::vector<int> map(nfaces + 1), lid(nfaces);

    // create maps from faces to DOFs
    int np(0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      int first = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      int shift(0);
      if (ndofs == 2) shift = UniqueIndexFaceToCells(*mesh_, f, c); 

      map[n] = np + shift;
      lid[n] = first + shift;
      np += ndofs;
    }
    map[nfaces] = np;
    
    // resize element matrices and right-hand side
    if (np != nfaces) {
      auto& Acell = local_op_->matrices[c];

      WhetStone::DenseMatrix Anew(np + 1, np + 1);
      Anew.PutScalar(0.0);

      for (int i = 0; i < nfaces + 1; ++i) {
        for (int j = 0; j < nfaces + 1; ++j) {
          Anew(map[i], map[j]) = Acell(i, j);
        }
      }

      local_op_->matrices[c] = Anew;
    }

    // shift gravity, k b g \grad z
    if (gravity_) {
      double zc = (mesh_->cell_centroid(c))[dim - 1];
      WhetStone::DenseMatrix& Wff = Wff_cells_[c];
      WhetStone::DenseVector v(nfaces), av(nfaces);

      double factor = rho_ * norm(g_);
      if (k_cell.get()) factor *= (*k_cell)[0][c];

      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        double zf = (mesh_->face_centroid(f))[dim - 1];
        v(n) = -(zf - zc) * factor;
      }

      Wff.Multiply(v, av, false);

      for (int n = 0; n < nfaces; n++) {
        rhs_face[0][lid[n]] += av(n); 
        rhs_cell[0][c] -= av(n); 
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: eliminate them or not. Finally we may add the essential BC
* the the system of equations as the trivial equations.
*
* Supported BCs on faces with many DOFs:
*   [Dirichlet]                         u = u0 
*   [Neumann]            -K(u) grad u . n = g0
*   [Mixed]        -K(u) grad u . n - c u = g1
****************************************************************** */
void PDE_DiffusionFracturedMatrix::ApplyBCs(
    bool primary, bool eliminate, bool essential_eqn)
{
  // apply diffusion type BCs to FACE-CELL system
  const std::vector<int>& bc_model_trial = bcs_trial_[0]->bc_model();
  const std::vector<int>& bc_model_test = bcs_test_[0]->bc_model();

  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  const std::vector<double>& bc_mixed = bcs_trial_[0]->bc_mixed();

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

  const auto& fmap = *cvs_->Map("face", true);

  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();
    
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    int nrows = Acell.NumRows();

    // un-roll multiple DOFs in a linear array
    std::vector<int> lid(nrows), mydof(nrows, 0), bctrial(nrows), bctest(nrows);
    std::vector<double> bcval(nrows), bcmix(nrows), bcarea(nrows);

    int np(0);
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int first = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      int shift(0);
      if (ndofs == 2) shift = UniqueIndexFaceToCells(*mesh_, f, c); 

      for (int k = 0; k < ndofs; ++k) {
        lid[np] = first + k;
        bctrial[np] = bc_model_trial[f];
        bctest[np] = bc_model_test[f];

        bcval[np] = bc_value[f];
        if (bc_mixed.size() > 0) bcmix[np] = bc_mixed[f];

        bcarea[np] = mesh_->face_area(f);
        if (k == shift) mydof[np] = 1;
        np++;
      }
    }

    // essential conditions for test functions
    for (int n = 0; n != np; ++n) {
      if (bctest[n] == OPERATOR_BC_DIRICHLET) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        for (int m = 0; m < nrows; m++) Acell(n, m) = 0.0;
      }
    }

    // conditions for trial functions
    for (int n = 0; n != np; ++n) {
      double value = bcval[n]; 
      if (bctrial[n] == OPERATOR_BC_DIRICHLET) {
        // make a copy of elemental matrix for post-processing
        if (flag) {
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }

        if (eliminate) { 
          for (int m = 0; m < np; m++) {
            rhs_face[0][lid[m]] -= Acell(m, n) * value;
            Acell(m, n) = 0.0;
          }

          rhs_cell[0][c] -= Acell(nrows - 1, n) * value;
          Acell(nrows - 1, n) = 0.0;
        }

        if (essential_eqn) {
          rhs_face[0][lid[n]] = value;
          Acell(n, n) = 1.0;
        }

      } else if (bctrial[n] == OPERATOR_BC_NEUMANN && primary && mydof[n] == 1) {
        rhs_face[0][lid[n]] -= value * bcarea[n];

      } else if (bctrial[n] == OPERATOR_BC_TOTAL_FLUX && primary && mydof[n] == 1) {
        rhs_face[0][lid[n]] -= value * bcarea[n];

      } else if (bctrial[n] == OPERATOR_BC_MIXED && primary && mydof[n] == 1) {
        if (flag) {  // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        rhs_face[0][lid[n]] -= value * bcarea[n];
        Acell(n, n) += bcmix[n] * bcarea[n];
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* WARNING: Since diffusive flux may be discontinuous (e.g. for
* Richards equation), we derive it only once (using flag) and in 
* exactly the same manner as other routines.
* **************************************************************** */
void PDE_DiffusionFracturedMatrix::UpdateFlux(
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<CompositeVector>& flux)
{
  Teuchos::RCP<const Epetra_MultiVector> k_cell = Teuchos::null;
  if (k_.get() && k_->HasComponent("cell")) k_cell = k_->ViewComponent("cell");

  // Initialize intensity in ghost faces.
  flux->PutScalarMasterAndGhosted(0.0);
  u->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& u_cell = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_face = *u->ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  int dim = mesh_->space_dimension();
  const auto& fmap = *cvs_->Map("face", true);

  int ndofs_owned = flux->ViewComponent("face")->MyLength();
  int ndofs_wghost = flux_data.MyLength();

  std::vector<int> hits(ndofs_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh_->cell_get_faces(c);
    const auto& dirs = mesh_->cell_get_face_dirs(c);
    int nfaces = faces.size();
    double zc = mesh_->cell_centroid(c)[dim - 1];

    // un-roll multiple DOFs in a linear array
    int nrows = 2 * nfaces;  // pessimistic estimate
    std::vector<int> lid(nrows), mydir(nrows), nohit(nrows, 0), map(nfaces);
    WhetStone::DenseVector v(nrows), av(nrows);

    int np(0);
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int first = fmap.FirstPointInElement(f);
      int ndofs = fmap.ElementSize(f);

      int shift(0);
      if (ndofs == 2) shift = UniqueIndexFaceToCells(*mesh_, f, c); 
      map[n] = np + shift;

      for (int k = 0; k < ndofs; ++k) {
        lid[np] = first + k;
        mydir[np] = dirs[n];
        if (ndofs > 1) nohit[np] = 1;

        v(np) = u_face[0][first + k];
        np++;
      }
    }
    v(np) = u_cell[0][c];

    v.Reshape(np + 1);
    av.Reshape(np + 1);

    if (local_op_->matrices_shadow[c].NumRows() == 0) { 
      local_op_->matrices[c].Multiply(v, av, false);
    } else {
      local_op_->matrices_shadow[c].Multiply(v, av, false);
    }

    if (gravity_) {
      double factor = rho_ * norm(g_);
      if (k_cell.get()) factor *= (*k_cell)[0][c];

      WhetStone::DenseVector w(nfaces), aw(nfaces);
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        double zf = (mesh_->face_centroid(f))[dim - 1];
        w(n) = -(zf - zc) * factor;
      }

      WhetStone::DenseMatrix& Wff = Wff_cells_[c];
      Wff.Multiply(w, aw, false);

      for (int n = 0; n < nfaces; n++) {
        av(map[n]) -= aw(n);
      }
    }

    // points of the master/slave interface require special logic
    for (int n = 0; n < np; n++) {
      int g = lid[n];
      if (g < ndofs_owned || nohit[n] == 1) {
        flux_data[0][g] -= av(n) * mydir[n];
        hits[g]++;
      }
      if (nohit[n] == 1) hits[g] = 1;
    }
  }

  for (int g = 0; g != ndofs_owned; ++g) {
    flux_data[0][g] /= hits[g];
  }

  flux->GatherGhostedToMaster();
}

}  // namespace Operators
}  // namespace Amanzi

