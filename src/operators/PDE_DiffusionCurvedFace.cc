/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "MatrixFE.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Diffusion_CurvedFace.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Schema.hh"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"
#include "Operator_Factory.hh"
#include "Operator_Schema.hh"
#include "UniqueLocalIndex.hh"

#include "PDE_DiffusionCurvedFace.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void
PDE_DiffusionCurvedFace::SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K)
{
  K_ = K;
  mass_matrices_initialized_ = false;
}


/* ******************************************************************
* Initialization of the operator: nonlinear coefficient.
****************************************************************** */
void
PDE_DiffusionCurvedFace::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                              const Teuchos::RCP<const CompositeVector>& dkdp)
{
  k_ = k;
  dkdp_ = dkdp;
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void
PDE_DiffusionCurvedFace::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                        const Teuchos::Ptr<const CompositeVector>& u)
{
  WhetStone::DenseMatrix Acell;
  WhetStone::MFD3D_Diffusion_CurvedFace mfd(plist_, mesh_);
  mfd.set_generalized_centroids(bf_);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  if (const_K_.rank() > 0) Kc = const_K_;

  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) Kc = (*K_)[c];
    mfd.StiffnessMatrix(c, Kc, Acell);
    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: eliminate them or not. Finally we may add the essntial BC
* the the system of equations as the trivial equations.
*
* NOTE 1. Nodal scheme handles only the case trialBC = testBC.
* NOTE 2. Jacobian term handles only trial BCs.
****************************************************************** */
void
PDE_DiffusionCurvedFace::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model_trial = bcs_trial_[0]->bc_model();
  const std::vector<int>& bc_model_test = bcs_test_[0]->bc_model();

  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();
  const std::vector<double>& bc_mixed = bcs_trial_[0]->bc_mixed();

  AMANZI_ASSERT(bc_model_trial.size() == nfaces_wghost);
  AMANZI_ASSERT(bc_value.size() == nfaces_wghost);

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_face = *global_op_->rhs()->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell");

  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    // essential conditions for test functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) {
        if (flag) { // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        for (int m = 0; m < nfaces + 1; m++) Acell(n, m) = 0.0;
      }
    }

    // conditions for trial functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      double value = bc_value[f];

      if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) {
        // make a copy of elemental matrix for post-processing
        if (flag) {
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }

        if (eliminate) {
          for (int m = 0; m < nfaces; m++) {
            rhs_face[0][faces[m]] -= Acell(m, n) * value;
            Acell(m, n) = 0.0;
          }

          rhs_cell[0][c] -= Acell(nfaces, n) * value;
          Acell(nfaces, n) = 0.0;
        }

        if (essential_eqn) {
          rhs_face[0][f] = value;
          Acell(n, n) = 1.0;
        }

      } else if (bc_model_trial[f] == OPERATOR_BC_NEUMANN && primary) {
        rhs_face[0][f] -= value * mesh_->face_area(f);

      } else if (bc_model_trial[f] == OPERATOR_BC_TOTAL_FLUX && primary) {
        rhs_face[0][f] -= value * mesh_->face_area(f);

      } else if (bc_model_trial[f] == OPERATOR_BC_MIXED && primary) {
        if (flag) { // make a copy of elemental matrix
          local_op_->matrices_shadow[c] = Acell;
          flag = false;
        }
        double area = mesh_->face_area(f);
        rhs_face[0][f] -= value * area;
        Acell(n, n) += bc_mixed[f] * area;
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Given pressures, reduce the problem to Lagrange multipliers.
****************************************************************** */
void
PDE_DiffusionCurvedFace::ModifyMatrices(const CompositeVector& u)
{
  if (local_op_schema_ !=
      (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
    std::cout << "Schema " << global_op_schema_ << " is not supported" << std::endl;
    AMANZI_ASSERT(0);
  }

  // populate the matrix
  const Epetra_MultiVector& u_c = *u.ViewComponent("cell");

  global_op_->rhs()->PutScalarGhosted(0.0);

  Epetra_MultiVector& rhs_f = *global_op_->rhs()->ViewComponent("face", true);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->cell_get_faces(c);
    int nfaces = faces.size();

    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      rhs_f[0][f] -= Acell(n, nfaces) * u_c[0][c];
      Acell(n, nfaces) = 0.0;
      Acell(nfaces, n) = 0.0;
    }
  }

  // Assemble all right-hand sides
  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* WARNING: Since diffusive flux may be discontinuous (e.g. for
* Richards equation), we derive it in exactly the same manner as
* in gravity routines.
* **************************************************************** */
void
PDE_DiffusionCurvedFace::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                                    const Teuchos::Ptr<CompositeVector>& flux)
{
  // Initialize intensity in ghost faces.
  flux->PutScalar(0.0);
  u->ScatterMasterToGhosted("face");

  const Epetra_MultiVector& u_cell = *u->ViewComponent("cell");
  const Epetra_MultiVector& u_face = *u->ViewComponent("face", true);
  Epetra_MultiVector& flux_data = *flux->ViewComponent("face", true);

  std::vector<int> hits(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& faces = mesh_->cell_get_faces(c);
    const auto& dirs = mesh_->cell_get_face_dirs(c);
    int nfaces = faces.size();

    WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
    for (int n = 0; n < nfaces; n++) { v(n) = u_face[0][faces[n]]; }
    v(nfaces) = u_cell[0][c];

    if (local_op_->matrices_shadow[c].NumRows() == 0) {
      local_op_->matrices[c].Multiply(v, av, false);
    } else {
      local_op_->matrices_shadow[c].Multiply(v, av, false);
    }

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (f < nfaces_owned) {
        flux_data[0][f] -= av(n) * dirs[n];
        hits[f]++;
      }
    }
  }

  for (int f = 0; f != nfaces_owned; ++f) { flux_data[0][f] /= hits[f]; }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.
****************************************************************** */
void
PDE_DiffusionCurvedFace::CreateMassMatrices_()
{
  WhetStone::MFD3D_Diffusion_CurvedFace mfd(plist_, mesh_);
  mfd.set_generalized_centroids(bf_);
  
  mfd.ModifyStabilityScalingFactor(factor_);

  WhetStone::DenseMatrix Wff;
  Wff_cells_.resize(ncells_owned);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  if (const_K_.rank() > 0) Kc = const_K_;

  for (int c = 0; c < ncells_owned; c++) {
    int ok;
    if (K_.get()) Kc = (*K_)[c];

    ok = mfd.MassMatrixInverse(c, Kc, Wff);
    Wff_cells_[c] = Wff;

    if (ok == 1) {
      Errors::Message msg("PDE_DiffusionCurvedFace: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }

  mass_matrices_initialized_ = true;
}


/* ******************************************************************
* Scale elemental inverse mass matrices. Use case is saturated flow.
****************************************************************** */
void
PDE_DiffusionCurvedFace::ScaleMassMatrices(double s)
{
  for (int c = 0; c < ncells_owned; c++) { Wff_cells_[c] *= s; }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void
PDE_DiffusionCurvedFace::Init_(Teuchos::ParameterList& plist)
{
  int d = mesh_->space_dimension();
  AMANZI_ASSERT(nfaces_owned > d * ncells_owned);

  if (global_op_ == Teuchos::null) { // create operator
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE;

    Schema schema(global_op_schema_);

    Operator_Factory factory;
    factory.set_mesh(mesh_);
    factory.set_plist(Teuchos::rcpFromRef(plist));
    factory.set_schema(schema);

    global_op_ = factory.CreateFromSchema();

  } else { // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  local_op_schema_ = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL;
  std::string name = "DiffusionCurvedFace: CELL_FACE+CELL";
  local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
  global_op_->OpPushBack(local_op_);

  // miscalleneous variables
  mass_matrices_initialized_ = false;
  K_ = Teuchos::null;
  k_ = Teuchos::null;
  dkdp_ = Teuchos::null;
  little_k_ = OPERATOR_LITTLE_K_NONE;

  // generalized face centroids are solution of a least-square problem
  // -- variable
  Teuchos::ParameterList tmp, tmp2;
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, d * d);

  Schema schema(AmanziMesh::CELL, d * d);
  Operator_Schema global_op(cvs, cvs, tmp, schema, schema);

  auto op = Teuchos::rcp(new Op_Face_Schema(schema, schema, mesh_));
  global_op.OpPushBack(op);

  // -- problem
  auto& rhs = *global_op.rhs();
  LSProblemSetup_(op->matrices, rhs);

  // -- solution
  global_op.SymbolicAssembleMatrix();
  global_op.AssembleMatrix();

  tmp.sublist("Hypre AMG")
     .set<std::string>("preconditioning method", "boomer amg").sublist("boomer amg parameters")
     .set<double>("tolerance", 0.0);
  tmp2.sublist("PCG").set<std::string>("iterative method", "pcg").sublist("pcg parameters")
     .set<Teuchos::Array<std::string>>("convergence criteria", std::vector<std::string>({"absolute residual", "relative residual", "make one iteration"}))
     .set<double>("error tolerance", 1e-12)
     .set<int>("maximum number of iterations", 400).sublist("verbose object")
     .set<std::string>("verbosity level", "medium");
  global_op.set_inverse_parameters("Hypre AMG", tmp, "PCG", tmp2);

  auto sol(rhs);
  global_op.ApplyInverse(rhs, sol);

  LSProblemPrimarySolution_(sol);

  // error analysis
  double err(0.0);
  for (int f = 0; f < nfaces_owned; ++f) {
    err += norm((*bf_)[f] - mesh_->face_centroid(f));
  }
  std::cout << "new face centroids: deviation=" << err / nfaces_owned << "\n\n";
}


/* ******************************************************************
* Setup problem A A^T w = I - A xref
****************************************************************** */
void
PDE_DiffusionCurvedFace::LSProblemSetup_(std::vector<WhetStone::DenseMatrix>& matrices,
                                         CompositeVector& rhs)
{
  int d = mesh_->space_dimension();
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  // matrix A A^T
  matrices.resize(nfaces_owned);

  for (int f = 0; f < nfaces_owned; ++f) {
    const auto& normal = mesh_->face_normal(f);
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    int nrows = ncells * d * d;
    WhetStone::DenseMatrix Aface(nrows, nrows);
    Aface.PutScalar(0.0);

    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) {
        double val = normal[i] * normal[j];
        for (int k = 0; k < d; ++k) {
          int i1 = k * d + i;
          int j1 = k * d + j;
          Aface(i1, j1) = val;
          if (ncells == 2) {
            int i2 = i1 + d * d;
            int j2 = j1 + d * d;
            Aface(i2, j2) = val;
            Aface(i1, j2) = -val;
            Aface(i2, j1) = -val;
          }
        } 
      }
    }

    matrices[f] = Aface;
  }

  // right-hand side is given by scaled identity matrices
  auto& rhs_c = *rhs.ViewComponent("cell", false);
  for (int c = 0; c < ncells_owned; ++c) {
    double vol = mesh_->cell_volume(c);
    int m(0);
    for (int i = 0; i < d; ++i) {
      for (int j = 0; j < d; ++j) rhs_c[m++][c] = ((i == j) ? vol : 0.0);
    }
  }

  // shift by vector A {xf} 
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      const auto& normal = mesh_->face_normal(f);
      const auto& xf = mesh_->face_centroid(f);

      int m(0);
      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) rhs_c[m++][c] -= normal[j] * xf[i] * dirs[n];
      }
    }
  }
}


/* ******************************************************************
* Recover solution vector x = A^T w
****************************************************************** */
void
PDE_DiffusionCurvedFace::LSProblemPrimarySolution_(const CompositeVector& sol)
{
  const auto& sol_c = *sol.ViewComponent("cell", true);
  sol.ScatterMasterToGhosted();

  // initialize generalized face centroids
  int d = mesh_->space_dimension();
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, d);
  auto bf = Teuchos::rcp(new CompositeVector(*cvs));

  auto bf_f = *bf->ViewComponent("face", true);
  for (int f = 0; f < nfaces_wghost; ++f) {
    const auto& xf = mesh_->face_centroid(f);
    for (int i = 0; i < d; ++i) bf_f[i][f] = xf[i];
  }

  // add correction to generalized face centroid
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      const auto& normal = mesh_->face_normal(f);

      int m(0);
      for (int i = 0; i < d; ++i) {
        for (int j = 0; j < d; ++j) bf_f[i][f] += sol_c[m++][c] * normal[j] * dirs[n];
      }
    }
  }

  bf->ScatterMasterToGhosted();

  // save to vector whocu could be shared with WhetStone
  bf_ = std::make_shared<std::vector<AmanziGeometry::Point>>(nfaces_wghost, AmanziGeometry::Point(d));
  for (int f = 0; f < nfaces_wghost; ++f) {
    for (int i = 0; i < d; ++i) (*bf_)[f][i] = bf_f[i][f];
  }
}

} // namespace Operators
} // namespace Amanzi
