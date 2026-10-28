/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

// Amanzi
#include "errors.hh"
#include "Inverse.hh"
#include "MatrixFE.hh"
#include "MFD3D_Diffusion.hh"
#include "InverseFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Face_Cell.hh"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"
#include "Operator_FaceCellSff.hh"

#include "PDE_DiffusionMFD.hh"

namespace Amanzi {
namespace Operators {


/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void
PDE_DiffusionMFD::SetTensorCoefficient(const Teuchos::RCP<const TensorVector>& K)
{
  transmissibility_initialized_ = false;
  K_ = K;
  if ((local_op_schema_ & OPERATOR_SCHEMA_BASE_CELL) && (K_ != Teuchos::null)) {
    AMANZI_ASSERT(K_->size() == ncells_owned);
  }
}

/* ******************************************************************
* Initialization of the operator: nonlinear coefficient.
****************************************************************** */
void
PDE_DiffusionMFD::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                       const Teuchos::RCP<const CompositeVector>& dkdp)
{
  if (k == Teuchos::null) {
    little_k_type_ = OPERATOR_LITTLE_K_NONE;
  } else {
    little_k_type_ = ParseLittleKType(plist_);
    AMANZI_ASSERT(little_k_type_ != OPERATOR_LITTLE_K_NONE);

    // compatibility checks
    if (little_k_type_ != OPERATOR_LITTLE_K_UPWIND) AMANZI_ASSERT(k->hasComponent("cell"));
    if (little_k_type_ != OPERATOR_LITTLE_K_STANDARD) AMANZI_ASSERT(k->hasComponent("face"));
    if (little_k_type_ == OPERATOR_LITTLE_K_DIVK_TWIN) AMANZI_ASSERT(k->hasComponent("twin"));
  }

  k_ = k;
  dkdp_ = dkdp;
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void
PDE_DiffusionMFD::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                 const Teuchos::Ptr<const CompositeVector>& u)
{
  // verify that mass matrices were initialized.
  if (!mass_matrices_initialized_) CreateMassMatrices_();

  bcs_applied_ = false;
  if (k_ != Teuchos::null) k_->scatterMasterToGhosted();

  if (!exclude_primary_terms_) {
    if (little_k_type_ == OPERATOR_LITTLE_K_NONE) {
      UpdateMatricesMixed_();
    } else {
      UpdateMatricesMixed_little_k_();
    }
  }
}


/* ******************************************************************
* Add stable approximation of Jacobian. It is done typically for
* the preconditioner.
****************************************************************** */
void
PDE_DiffusionMFD::UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                                 const Teuchos::Ptr<const CompositeVector>& u,
                                                 double scalar_factor)
{
  //  add Newton-type corrections
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE) {
    if (dkdp_ != Teuchos::null) dkdp_->scatterMasterToGhosted();
    AddNewtonCorrectionCell_(flux, u, scalar_factor);
  }
}


/* ******************************************************************
* Basic routine for each operator: creation of elemental matrices.
****************************************************************** */
void
PDE_DiffusionMFD::UpdateMatricesMixed_()
{
  DenseMatrix_Vector& A = local_op_->A;
  auto& Wff_cells = Wff_cells_;

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::UpdateMatricesMixed_", ncells_owned, KOKKOS_LAMBDA(const int& c) {
      auto Wff = Wff_cells[c];
      int nfaces = Wff.NumRows();

      auto Acell = A[c];

      // create stiffness matrix by elimination of the mass matrix
      double matsum = 0.0;
      for (int n = 0; n < nfaces; n++) {
        double rowsum = 0.0;
        for (int m = 0; m < nfaces; m++) {
          double tmp = Wff(n, m);
          rowsum += tmp;
          Acell(n, m) = tmp;
        }
        Acell(n, nfaces) = -rowsum;
        Acell(nfaces, n) = -rowsum;
        matsum += rowsum;
      }
      Acell(nfaces, nfaces) = matsum;
    });
}


/* ******************************************************************
* Helper function for staging kr
****************************************************************** */
void
PDE_DiffusionMFD::Preallocate_little_k_(bool init)
{
  // pre-stage kr face and cell quantities into work space
  if (kr_cells_.size() != Wff_cells_.size()) {
    kr_cells_ = DenseVector_Vector(Wff_cells_.size());
    for (int i = 0; i != Wff_cells_.size(); ++i) {
      kr_cells_.set_shape(i, Wff_cells_.at_host(i).NumRows());
    }
    kr_cells_.Init();
    if (init) kr_cells_.putScalar(1.0);
  }
}

/* ******************************************************************
* Basic routine for each operator: creation of elemental matrices.
****************************************************************** */
void
PDE_DiffusionMFD::UpdateMatricesMixed_little_k_()
{
  if (k_ != Teuchos::null) k_->scatterMasterToGhosted();
  CompositeVector::cView_type k_cell, k_face, k_twin;
  if (k_ != Teuchos::null) {
    if (k_->hasComponent("cell")) k_cell = k_->viewComponent("cell", false);
    if (k_->hasComponent("face")) k_face = k_->viewComponent("face", true);
    if (k_->hasComponent("twin")) k_twin = k_->viewComponent("twin", true);
  }

  DenseMatrix_Vector& A = local_op_->A;
  const AmanziMesh::MeshCache& m = mesh_->getCache();
  Preallocate_little_k_();

  auto& Wff_cells = Wff_cells_;
  auto& kr_cells = kr_cells_;
  int little_k_type(little_k_type_);

  Kokkos::parallel_for(
    "PDE_DiffusionMFD::UpdateMatricesMixed_little_k_", ncells_owned, KOKKOS_LAMBDA(const int& c) {
      auto faces = m.getCellFaces(c);
      int nfaces = faces.extent(0);

      // set up kr
      auto kr = kr_cells[c];
      kr.putScalar(1.0);

      // -- chefs recommendation: SPD discretization with upwind
      if (little_k_type == OPERATOR_LITTLE_K_DIVK) {
        if (k_face.extent(0) > 0)
          for (int n = 0; n < nfaces; n++) kr(n) = k_face(faces[n], 0);
        if (k_cell.extent(0) > 0) kr(nfaces) = k_cell(c, 0);

        // -- new scheme: SPD discretization with upwind and equal splitting
      } else if (little_k_type == OPERATOR_LITTLE_K_DIVK_BASE) {
        if (k_face.extent(0) > 0)
          for (int n = 0; n < nfaces; n++) kr(n) = sqrt(k_face(faces[n], 0));

        // -- same as above but remains second-order for dicontinuous coefficients
      } else if (little_k_type == OPERATOR_LITTLE_K_DIVK_TWIN) {
        if (k_face.extent(0) > 0) {
          for (int n = 0; n < nfaces; n++) {
            int f = faces[n];
            auto cells = m.getFaceCells(f);
            kr(n) = (c == cells[0]) ? k_face(f, 0) : k_twin(f, 0);
          }
        }
        if (k_cell.extent(0) > 0) kr(nfaces) = k_cell(c, 0);

        // -- the second most popular choice: classical upwind
      } else if (little_k_type == OPERATOR_LITTLE_K_UPWIND) {
        if (k_face.extent(0) > 0)
          for (int n = 0; n < nfaces; n++) kr(n) = k_face(faces[n], 0);

      } else if (little_k_type == OPERATOR_LITTLE_K_STANDARD) {
        if (k_cell.extent(0) > 0)
          for (int n = 0; n < nfaces; n++) kr(n) = k_cell(c, 0);
      }

      auto Wff = Wff_cells[c];
      auto Acell = A[c];

      if (little_k_type & OPERATOR_LITTLE_K_DIVK_BASE) {
        double matsum = 0.0;
        for (int n = 0; n < nfaces; n++) {
          double rowsum = 0.0;
          for (int m = 0; m < nfaces; m++) {
            double tmp = Wff(n, m) * kr(n) * kr(m) / kr(nfaces);
            rowsum += tmp;
            Acell(n, m) = tmp;
          }

          Acell(n, nfaces) = -rowsum;
          Acell(nfaces, n) = -rowsum;
          matsum += rowsum;
        }

        Acell(nfaces, nfaces) = matsum;

      } else {
        double matsum = 0.0;
        for (int n = 0; n < nfaces; n++) {
          double rowsum = 0.0;
          for (int m = 0; m < nfaces; m++) {
            double tmp = Wff(n, m) * kr(n);
            rowsum += tmp;
            Acell(n, m) = tmp;
          }

          Acell(n, nfaces) = -rowsum;
          matsum += rowsum;
        }
        Acell(nfaces, nfaces) = matsum;

        for (int n = 0; n < nfaces; n++) {
          double colsum = 0.0;
          for (int m = 0; m < nfaces; m++) colsum += Acell(m, n);
          Acell(nfaces, n) = -colsum;
        }
      }
    });
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
PDE_DiffusionMFD::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  bcs_applied_ = true;

  if (!exclude_primary_terms_) {
    AMANZI_ASSERT(bcs_trial_.size() == 1);
    AMANZI_ASSERT(bcs_test_.size() == 1);
    ApplyBCs_Mixed_(bcs_trial_[0].ptr(), bcs_test_[0].ptr(), primary, eliminate, essential_eqn);
  }
  if (jac_op_.get()) ApplyBCsJacobian();
}

void
PDE_DiffusionMFD::ApplyBCsJacobian()
{
  const auto bc_model = bcs_trial_[0]->bc_model();
  DenseMatrix_Vector& A = jac_op_->A;
  Kokkos::parallel_for(
    "PDE_DiffusionMFD::ApplyBCs on Jacobian", nfaces_owned, KOKKOS_LAMBDA(const int& f) {
      if (bc_model[f] == OPERATOR_BC_NEUMANN || bc_model[f] == OPERATOR_BC_TOTAL_FLUX) {
        auto Aface = A[f];
        Aface *= 0.0;
      }
    });
}


/* ******************************************************************
* Apply BCs on face values.
****************************************************************** */
void
PDE_DiffusionMFD::ApplyBCs_Mixed_(const Teuchos::Ptr<const BCs>& bc_trial,
                                  const Teuchos::Ptr<const BCs>& bc_test,
                                  bool primary,
                                  bool eliminate,
                                  bool essential_eqn)
{
  global_op_->rhs()->putScalarGhosted(0.0);

  { // context for views
    // apply diffusion type BCs to FACE-CELL system
    const AmanziMesh::MeshCache& m = mesh_->getCache();
    DenseMatrix_Vector& A = local_op_->A;

    const auto bc_model_trial = bc_trial->bc_model();
    const auto bc_model_test = bc_test->bc_model();

    const auto bc_value = bc_trial->bc_value();
    const auto bc_mixed = bc_trial->bc_mixed();

    AMANZI_ASSERT(bc_model_trial.size() == nfaces_wghost);
    AMANZI_ASSERT(bc_value.size() == nfaces_wghost);

    const auto rhs_face = global_op_->rhs()->viewComponent("face", true);
    const auto rhs_cell = global_op_->rhs()->viewComponent("cell");

    Kokkos::parallel_for(
      "PDE_DiffusionMFD::ApplyBCs_Mixed_", ncells_owned, KOKKOS_LAMBDA(const int& c) {
        auto faces = m.getCellFaces(c);
        int nfaces = faces.size();

        auto Acell = A[c];
        // essential conditions for test functions
        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          if (bc_model_test[f] == OPERATOR_BC_DIRICHLET) {
            for (int m = 0; m < nfaces + 1; m++) Acell(n, m) = 0.0;
          }
        }

        // conditions for trial functions
        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          double value = bc_value[f];

          if (bc_model_trial[f] == OPERATOR_BC_DIRICHLET) {
            if (eliminate) {
              for (int m = 0; m < nfaces; m++) {
                Kokkos::atomic_add(&rhs_face(faces[m], 0), -Acell(m, n) * value);
                Acell(m, n) = 0.0;
              }

              rhs_cell(c, 0) -= Acell(nfaces, n) * value;
              Acell(nfaces, n) = 0.0;
            }

            if (essential_eqn) {
              rhs_face(f, 0) = value;
              Acell(n, n) = 1.0;
            }
          } else if (bc_model_trial[f] == OPERATOR_BC_NEUMANN && primary) {
            // need not be atomic if f is boundary face?
            Kokkos::atomic_add(&rhs_face(f, 0), -value * m.getFaceArea(f));
          } else if (bc_model_trial[f] == OPERATOR_BC_TOTAL_FLUX && primary) {
            Kokkos::atomic_add(&rhs_face(f, 0), -value * m.getFaceArea(f));
          } else if (bc_model_trial[f] == OPERATOR_BC_MIXED && primary) {
            double area = m.getFaceArea(f);
            Kokkos::atomic_add(&rhs_face(f, 0), -value * area);
            Acell(n, n) += bc_mixed[f] * area;
          }
        }
      });
  }
  global_op_->rhs()->gatherGhostedToMaster("face", Tpetra::ADD);
}


/* ******************************************************************
* Modify operator by adding upwind approximation of Newton corection.
* A special care should be taken later to deal with the case
* where kf < 0, i.e. for energy when: div (qh) = div (h k grad p),
* where h is enthalpy and can be negative. I think that the current
* treatment is inadequate.
****************************************************************** */
void
PDE_DiffusionMFD::AddNewtonCorrectionCell_(const Teuchos::Ptr<const CompositeVector>& flux,
                                           const Teuchos::Ptr<const CompositeVector>& u,
                                           double scalar_factor)
{
  // hack: ignore correction if no flux provided.
  if (flux == Teuchos::null) return;

  // Correction is zero for linear problems
  if (k_ == Teuchos::null || dkdp_ == Teuchos::null) return;

  // only works on upwinded methods
  if (little_k_type_ == OPERATOR_UPWIND_NONE) return;

  { // context for views
    const AmanziMesh::MeshCache& m = mesh_->getCache();
    DenseMatrix_Vector& A = jac_op_->A;

    const auto dkdp_f = dkdp_->viewComponent("face");
    const auto flux_f = flux->viewComponent("face");
    const auto kf = k_->viewComponent("face");

    Kokkos::parallel_for(
      "PDE_DiffusionMFD::AddNewtonCorrectionCell_", nfaces_owned, KOKKOS_LAMBDA(const int& f) {
        // populate the local matrices
        auto cells = m.getFaceCells(f);
        int ncells = cells.size();

        auto Aface = A[f];
        Aface.putScalar(0.);

        // We use the upwind discretization of the generalized flux.
        double v = Kokkos::abs(kf(f, 0)) > 0.0 ? flux_f(f, 0) * dkdp_f(f, 0) / kf(f, 0) : 0.0;
        double vmod = Kokkos::abs(v);

        // prototype for future limiters (external or internal ?)
        vmod *= scalar_factor;

        // define the upwind cell, index i in this case
        int i, dir, c1;
        c1 = cells[0];
        const AmanziGeometry::Point normal = m.getFaceNormal(f, c1, &dir);
        i = (v * dir >= 0.0) ? 0 : 1;

        if (ncells == 2) {
          Aface(i, i) = vmod;
          Aface(1 - i, i) = -vmod;
        } else if (i == 0) {
          Aface(0, 0) = vmod;
        }
      });
  }
}


/* ******************************************************************
* WARNING: Since diffusive flux may be discontinuous (e.g. for
* Richards equation), we derive it in exactly the same manner as
* in gravity routines.
* **************************************************************** */
void
PDE_DiffusionMFD::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                             const Teuchos::Ptr<CompositeVector>& flux)
{
  if (bcs_applied_) {
    Errors::Message msg("PDE_DiffusionMFD::UpdateFlux: Developer Error -- UpdateFlux() cannot be "
                        "called after BCs have been applied.");
    Exceptions::amanzi_throw(msg);
  }

  flux->putScalarMasterAndGhosted(0.0);
  u->scatterMasterToGhosted("face");
  CompositeVector_<int> hits(flux->getMap(), InitMode::ZERO);

  if (k_ != Teuchos::null) {
    if (k_->hasComponent("face")) k_->scatterMasterToGhosted("face");
  }

  if (local_op_->v.size() != local_op_->A.size()) { local_op_->PreallocateWorkVectors(); }

  {
    const auto u_cell = u->viewComponent("cell");
    const auto u_face = u->viewComponent("face", true);

    auto flux_data = flux->viewComponent("face", true);
    auto hits_data = hits.viewComponent("face", true);

    const AmanziMesh::MeshCache& m = mesh_->getCache();

    auto local_A = local_op_->A;
    auto local_Av = local_op_->Av;
    auto local_v = local_op_->v;

    Kokkos::parallel_for(
      "PDE_DiffusionMFD::UpdateFlux", ncells_owned, KOKKOS_LAMBDA(const int& c) {
        auto [faces, dirs] = m.getCellFacesAndDirections(c);
        int nfaces = faces.size();

        auto lv = local_v[c];
        auto lAv = local_Av[c];
        auto lA = local_A[c];

        for (int n = 0; n < nfaces; n++) { lv(n) = u_face(faces[n], 0); }
        lv(nfaces) = u_cell(c, 0);

        lA.Multiply(lv, lAv, false);

        for (int n = 0; n < nfaces; n++) {
          int f = faces[n];
          Kokkos::atomic_add(&flux_data(f, 0), -lAv(n) * dirs[n]);
          Kokkos::atomic_inc(&hits_data(f, 0));
        }
      });
  }

  flux->gatherGhostedToMaster("face", Tpetra::ADD);
  hits.gatherGhostedToMaster("face", Tpetra::ADD);

  {
    auto flux_data = flux->viewComponent("face", false);
    auto hits_data = hits.viewComponent("face", false);
    Kokkos::parallel_for(
      "PDE_DiffusionMFD::UpdateFlux", nfaces_owned, KOKKOS_LAMBDA(const int& f) {
        flux_data(f, 0) /= hits_data(f, 0);
      });
  }
}


/* ******************************************************************
* Calculate elemental inverse mass matrices.
****************************************************************** */
void
PDE_DiffusionMFD::CreateMassMatrices_()
{
  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix<> Wff;

  // not supporting surface mesh
  assert(mesh_->getManifoldDimension() == mesh_->getSpaceDimension());

  std::vector<WhetStone::DenseMatrix<>> tmp_Wff_cells;
  tmp_Wff_cells.resize(ncells_owned);

  WhetStone::Tensor<> Kc(mesh_->getSpaceDimension(), 1);
  Kc(0, 0) = 1.0;

  K_->update_entries_host();

  for (int c = 0; c < ncells_owned; c++) {
    int ok(WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED);
    if (K_.get()) Kc = K_->at_host(c);

    // For problems with degenerate coefficients we should skip WhetStone.
    assert(Kc.Trace() != 0.0);

    int method = mfd_primary_;
    ok = WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED;

    // try primary and then secondary discretization methods.
    if (method == WhetStone::DIFFUSION_HEXAHEDRA_MONOTONE) {
      assert(false);
      ok = mfd.MassMatrixInverseMMatrixHex(c, Kc, Wff);
      method = mfd_secondary_;
    } else if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_MONOTONICITY) {
      assert(false);
      ok = mfd.MassMatrixInverseMMatrix(c, Kc, Wff);
      method = mfd_secondary_;
    }

    if (ok != WhetStone::WHETSTONE_ELEMENTAL_MATRIX_OK) {
      if (method == WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY) {
        ok = mfd.MassMatrixInverseOptimized(c, Kc, Wff);
      } else if (method == WhetStone::DIFFUSION_TPFA) {
        ok = mfd.MassMatrixInverseTPFA(c, Kc, Wff);
      } else if (method == WhetStone::DIFFUSION_SUPPORT_OPERATOR) {
        assert(false);
        ok = mfd.MassMatrixInverseSO(c, Kc, Wff);
      } else if (method == WhetStone::DIFFUSION_POLYHEDRA_SCALED) {
        if (K_symmetric_) {
          ok = mfd.MassMatrixInverse(c, Kc, Wff);
        } else {
          ok = mfd.MassMatrixInverseNonSymmetric(c, Kc, Wff);
        }
      }
    }

    tmp_Wff_cells[c].reshape(Wff.NumRows(), Wff.NumCols());
    tmp_Wff_cells[c].assign(Wff);

    if (ok == WhetStone::WHETSTONE_ELEMENTAL_MATRIX_FAILED) {
      Errors::Message msg("PDE_DiffusionMFD: unexpected failure in WhetStone.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // Copy data to GPU
  Wff_cells_.Init(tmp_Wff_cells);
  mass_matrices_initialized_ = true;
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void
PDE_DiffusionMFD::ParsePList_()
{
  // Determine discretization
  std::string primary = plist_.get<std::string>("discretization primary");
  std::string secondary = plist_.get<std::string>("discretization secondary", primary);
  K_symmetric_ = (plist_.get<std::string>("diffusion tensor", "symmetric") == "symmetric");

  // Primary discretization methods
  if (primary == "mfd: two-point flux approximation") {
    mfd_primary_ = WhetStone::DIFFUSION_TPFA;
  } else if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  } else if (primary == "mfd: optimized for sparsity") {
    mfd_primary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  } else if (primary == "mfd: support operator") {
    mfd_primary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: primary discretization method \"" << primary
        << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Secondary discretization methods
  if (secondary == "mfd: two-point flux approximation") {
    mfd_secondary_ = WhetStone::DIFFUSION_TPFA;
  } else if (secondary == "mfd: default") {
    mfd_secondary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  } else if (secondary == "mfd: optimized for sparsity") {
    mfd_secondary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  } else if (secondary == "mfd: support operator") {
    mfd_secondary_ = WhetStone::DIFFUSION_SUPPORT_OPERATOR;
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: secondary discretization method \"" << secondary
        << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Define stencil for the MFD diffusion method.
  std::vector<std::string> names;
  if (plist_.isParameter("schema")) {
    names = plist_.get<Teuchos::Array<std::string>>("schema").toVector();
  } else {
    names.resize(2);
    names[0] = "face";
    names[1] = "cell";
    plist_.set<Teuchos::Array<std::string>>("schema", names);
  }

  // map schema into DOFS
  local_op_schema_ = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL;
  schema_prec_dofs_ = OPERATOR_SCHEMA_DOFS_FACE | OPERATOR_SCHEMA_DOFS_CELL;
}


/* ******************************************************************
 * Virtual constructor, sets up global op, cvs.
 ****************************************************************** */
void
PDE_DiffusionMFD::Init()
{
  // create or check the existing Operator
  int global_op_schema = schema_prec_dofs_;
  if (global_op_ == Teuchos::null) {
    global_op_schema_ = global_op_schema;

    // build the CVS from the global schema
    auto cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

    global_op_ = Teuchos::rcp(new Operator_FaceCell(cvs->CreateSpace(), plist_));

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->getDomainMap()->getMesh();
  }

  // Do we need to exclude the primary terms for a Jacobian-only off-diagonal?
  exclude_primary_terms_ = plist_.get<bool>("exclude primary terms", false);

  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms_) {
    std::string name = "Diffusion: CELL_FACE+CELL";
    local_op_ = Teuchos::rcp(new Op_Cell_FaceCell(name, mesh_));
    global_op_->OpPushBack(local_op_);
  }

  // scalar coefficient
  little_k_type_ = ParseLittleKType(plist_);

  // Do we need to calculate Newton correction terms?
  std::string jacobian = plist_.get<string>("Newton correction", "none");
  if (jacobian == "none") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  } else if (jacobian == "approximate Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;

    // create a local op
    jac_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;
    std::string opname("Jacobian FACE_CELL");
    jac_op_ = Teuchos::rcp(new Op_Face_Cell(opname, mesh_));
    global_op_->OpPushBack(jac_op_);
  } else {
    Errors::Message msg;
    msg << "PDE_DiffusionMFD: unrecognized Jacobian option \"" << jacobian << "\".";
    Exceptions::amanzi_throw(msg);
  }
  
  // miscalleneous variables
  K_ = Teuchos::null;
  k_ = Teuchos::null;
  dkdp_ = Teuchos::null;
}


} // namespace Operators
} // namespace Amanzi
