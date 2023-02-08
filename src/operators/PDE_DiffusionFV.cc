/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// Operators
#include "Op.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"
#include "Op_Face_Cell.hh"
#include "OperatorDefs.hh"
#include "Operator_Cell.hh"
#include "PDE_DiffusionFV.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
****************************************************************** */
void
PDE_DiffusionFV::Init_(Teuchos::ParameterList& plist)
{
  // Define stencil for the FV diffusion method.
  local_op_schema_ = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_CELL;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    global_op_schema_ = OPERATOR_SCHEMA_DOFS_CELL;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist, global_op_schema_));

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // Do we need to exclude the primary terms?
  exclude_primary_terms_ = plist.get<bool>("exclude primary terms", false);

  // create the local Op and register it with the global Operator
  if (!exclude_primary_terms_) {
    if (plist.get<bool>("surface operator", false)) {
      std::string name = "Diffusion: FACE_CELL Surface";
      local_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
    } else {
      std::string name = "Diffusion: FACE_CELL";
      local_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
    }
    global_op_->OpPushBack(local_op_);
  }

  // upwind options
  Errors::Message msg;
  std::string uwname = plist.get<std::string>("nonlinear coefficient", "upwind: face");
  if (uwname == "none") {
    little_k_ = OPERATOR_LITTLE_K_NONE;

  } else if (uwname == "upwind: face") {
    little_k_ = OPERATOR_LITTLE_K_UPWIND;

  } else if (uwname == "divk: cell-face" || uwname == "divk: cell-grad-face-twin" ||
             uwname == "standard: cell") {
    msg << "DiffusionFV: \"" << uwname << "\" upwinding not supported.";
    Exceptions::amanzi_throw(msg);

  } else {
    msg << "DiffusionFV: unknown upwind scheme specified.";
    Exceptions::amanzi_throw(msg);
  }

  // Do we need to calculate Newton correction terms?
  newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;

  // DEPRECATED INPUT -- remove this error eventually --etc
  if (plist.isParameter("newton correction")) {
    msg << "DiffusionFV: DEPRECATED: \"newton correction\" has been removed in favor of \"Newton "
           "correction\"";
    Exceptions::amanzi_throw(msg);
  }

  std::string jacobian = plist.get<std::string>("Newton correction", "none");
  if (jacobian == "none") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_NONE;
  } else if (jacobian == "true Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_TRUE;
  } else if (jacobian == "approximate Jacobian") {
    newton_correction_ = OPERATOR_DIFFUSION_JACOBIAN_APPROXIMATE;
    msg << "DiffusionFV: \"approximate Jacobian\" not supported, use \"true Jacobian\".";
    Exceptions::amanzi_throw(msg);
  } else {
    msg << "DiffusionFV: invalid parameter \"" << jacobian
        << "\" for option \"Newton correction\" -- valid are: \"none\", \"true Jacobian\"";
    Exceptions::amanzi_throw(msg);
  }


  if (newton_correction_ != OPERATOR_DIFFUSION_JACOBIAN_NONE) {
    if (plist.get<bool>("surface operator", false)) {
      std::string name = "Diffusion: FACE_CELL Surface Jacobian terms";
      jac_op_ = Teuchos::rcp(new Op_SurfaceFace_SurfaceCell(name, mesh_));
    } else {
      std::string name = "Diffusion: FACE_CELL Jacobian terms";
      jac_op_ = Teuchos::rcp(new Op_Face_Cell(name, mesh_));
    }

    global_op_->OpPushBack(jac_op_);
  }

  // solution-independent data
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  transmissibility_ = Teuchos::rcp(new CompositeVector(cvs, true));
}


/* ******************************************************************
* Setup methods: scalar coefficients
****************************************************************** */
void
PDE_DiffusionFV::SetTensorCoefficient(const Teuchos::RCP<const std::vector<WhetStone::Tensor>>& K)
{
  transmissibility_initialized_ = false;
  K_ = K;
}


/* ******************************************************************
* Setup methods: krel and deriv -- must be called after calling a
* setup with K absolute
****************************************************************** */
void
PDE_DiffusionFV::SetScalarCoefficient(const Teuchos::RCP<const CompositeVector>& k,
                                      const Teuchos::RCP<const CompositeVector>& dkdp)
{
  transmissibility_initialized_ = false;
  k_ = k;
  dkdp_ = dkdp;

  if (k_ != Teuchos::null) {
    AMANZI_ASSERT(k_->HasComponent("face"));
    AMANZI_ASSERT(k_->Map().Ghosted());

    // NOTE: it seems that Amanzi passes in a cell based kr which is then
    // ignored, and assumed = 1.  This seems dangerous to me. --etc
    // AMANZI_ASSERT(!k_->HasComponent("cell"));
  }
  if (dkdp_ != Teuchos::null) {
    AMANZI_ASSERT(dkdp_->HasComponent("cell"));
    AMANZI_ASSERT(dkdp_->Map().Ghosted());
  }
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces.
****************************************************************** */
void
PDE_DiffusionFV::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                const Teuchos::Ptr<const CompositeVector>& u)
{
  if (!transmissibility_initialized_) ComputeTransmissibility_();

  if (!exclude_primary_terms_) {
    const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);

    // preparing upwind data
    Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
    if (k_ != Teuchos::null) {
      k_->ScatterMasterToGhosted("face");
      if (k_->HasComponent("face")) k_face = k_->ViewComponent("face", true);
    }

    // updating matrix blocks
    AmanziMesh::Entity_ID_View cells, faces;

    for (int f = 0; f != nfaces_owned; ++f) {
      cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int ncells = cells.size();

      WhetStone::DenseMatrix Aface(ncells, ncells);
      Aface = 0.0;

      double tij = trans_face[0][f] * (k_face.get() ? (*k_face)[0][f] : 1.0);
      for (int i = 0; i != ncells; ++i) {
        Aface(i, i) = tij;
        for (int j = i + 1; j != ncells; ++j) {
          Aface(i, j) = -tij;
          Aface(j, i) = -tij;
        }
      }
      local_op_->matrices[f] = Aface;
    }
  }
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces with Newton information
****************************************************************** */
void
PDE_DiffusionFV::UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                                const Teuchos::Ptr<const CompositeVector>& u,
                                                double scalar_factor)
{
  // Add derivatives to the matrix (Jacobian in this case)
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_TRUE && u.get()) {
    AMANZI_ASSERT(u != Teuchos::null);
    AnalyticJacobian_(*u);
  }
}

void
PDE_DiffusionFV::UpdateMatricesNewtonCorrection(const Teuchos::Ptr<const CompositeVector>& flux,
                                                const Teuchos::Ptr<const CompositeVector>& u,
                                                const Teuchos::Ptr<const CompositeVector>& factor)
{
  // Add derivatives to the matrix (Jacobian in this case)
  if (newton_correction_ == OPERATOR_DIFFUSION_JACOBIAN_TRUE && u.get()) {
    AMANZI_ASSERT(u != Teuchos::null);
    AnalyticJacobian_(*u);
  }
}


/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void
PDE_DiffusionFV::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);

  AMANZI_ASSERT(bcs_trial_.size() > 0);
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  Teuchos::RCP<const Epetra_MultiVector> k_face = Teuchos::null;
  if (k_ != Teuchos::null) { k_face = k_->ViewComponent("face", true); }

  if (!exclude_primary_terms_) {
    Epetra_MultiVector& rhs_cell = *global_op_->rhs()->ViewComponent("cell", true);

    AmanziMesh::Entity_ID_View cells;

    for (int f = 0; f < nfaces_owned; f++) {
      if (bc_model[f] != OPERATOR_BC_NONE) {
        cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        int c = cells[0];

        if (bc_model[f] == OPERATOR_BC_DIRICHLET && primary) {
          rhs_cell[0][c] += bc_value[f] * trans_face[0][f] * (k_face.get() ? (*k_face)[0][f] : 1.0);
        } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
          local_op_->matrices_shadow[f] = local_op_->matrices[f];
          local_op_->matrices[f](0, 0) = 0.0;

          if (primary) rhs_cell[0][c] -= bc_value[f] * mesh_->getFaceArea(f);
        }
      }
    }
  }

  if (jac_op_ != Teuchos::null) {
    AMANZI_ASSERT(bc_model.size() == nfaces_wghost);
    for (int f = 0; f != nfaces_owned; ++f) {
      WhetStone::DenseMatrix& Aface = jac_op_->matrices[f];

      if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        jac_op_->matrices_shadow[f] = Aface;
        Aface *= 0.0;
      }
    }
  }
}


/* ******************************************************************
* Calculate mass flux from cell-centered data
****************************************************************** */
void
PDE_DiffusionFV::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                            const Teuchos::Ptr<CompositeVector>& darcy_mass_flux)
{
  const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  solution->ScatterMasterToGhosted("cell");
  if (k_ != Teuchos::null) k_->ScatterMasterToGhosted("face");

  const Teuchos::Ptr<const Epetra_MultiVector> Krel_face =
    k_.get() ? k_->ViewComponent("face", false).ptr() : Teuchos::null;

  const Epetra_MultiVector& p = *solution->ViewComponent("cell", true);
  Epetra_MultiVector& flux = *darcy_mass_flux->ViewComponent("face", false);

  AmanziMesh::Entity_ID_View cells;

  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        double value = bc_value[f];
        flux[0][f] = dirs[n] * trans_face[0][f] * (p[0][c] - value);
        if (Krel_face.get()) flux[0][f] *= (*Krel_face)[0][f];

      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        double value = bc_value[f];
        double area = mesh_->getFaceArea(f);
        flux[0][f] = dirs[n] * value * area;

      } else {
        if (f < nfaces_owned && !flag[f]) {
          cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
          if (cells.size() <= 1) {
            Errors::Message msg("Flow PK: These boundary conditions are not supported by FV.");
            Exceptions::amanzi_throw(msg);
          }
          int c1 = cells[0];
          int c2 = cells[1];
          if (c == c1) {
            flux[0][f] = dirs[n] * trans_face[0][f] * (p[0][c1] - p[0][c2]);
          } else {
            flux[0][f] = dirs[n] * trans_face[0][f] * (p[0][c2] - p[0][c1]);
          }
          if (Krel_face.get()) flux[0][f] *= (*Krel_face)[0][f];
          flag[f] = 1;
        }
      }
    }
  }
}


/* ******************************************************************
* Scale face-based matrices.
****************************************************************** */
void
PDE_DiffusionFV::ScaleMatricesColumns(const CompositeVector& s)
{
  AmanziMesh::Entity_ID_List cells;

  const auto& s_c = *s.ViewComponent("cell");

  for (int f = 0; f < nfaces_owned; ++f) {
    WhetStone::DenseMatrix& Aface = local_op_->matrices[f];

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; ++n) {
      double factor = s_c[0][cells[n]];
      for (int m = 0; m < ncells; ++m) { Aface(m, n) *= factor; }
    }
  }
}


/* ******************************************************************
* Computation the part of the Jacobian which depends on derivatives
* of the relative permeability wrt to capillary pressure. They must
* be added to the existing matrix structure.
****************************************************************** */
void
PDE_DiffusionFV::AnalyticJacobian_(const CompositeVector& u)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<double>& bc_value = bcs_trial_[0]->bc_value();

  u.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);

  AmanziMesh::Entity_ID_View cells, faces;
  double dkdp[2], pres[2];

  dkdp_->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& dKdP_cell = *dkdp_->ViewComponent("cell", true);
  AMANZI_ASSERT(dKdP_cell.MyLength() == ncells_wghost);

  Teuchos::RCP<const Epetra_MultiVector> dKdP_face;
  if (dkdp_->HasComponent("face")) { dKdP_face = dkdp_->ViewComponent("face", true); }

  for (int f = 0; f != nfaces_owned; ++f) {
    cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int mcells = cells.size();

    WhetStone::DenseMatrix Aface(mcells, mcells);
    for (int n = 0; n < mcells; n++) {
      int c1 = cells[n];
      pres[n] = uc[0][c1];
      dkdp[n] = dKdP_cell[0][c1];
    }

    if (mcells == 1) dkdp[1] = dKdP_face.get() ? (*dKdP_face)[0][f] : 0.;

    // find the face direction from cell 0 to cell 1
    const auto& [cfaces,fdirs] = mesh_->getCellFacesAndDirections(cells[0]);

    int f_index = std::find(cfaces.begin(), cfaces.end(), f) - cfaces.begin();
    ComputeJacobianLocal_(mcells, f, fdirs[f_index], bc_model[f], bc_value[f], pres, dkdp, Aface);
    jac_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Computation of a local submatrix of the analytical Jacobian
* (its nonlinear part) on face f.
****************************************************************** */
void
PDE_DiffusionFV::ComputeJacobianLocal_(int mcells,
                                       int f,
                                       int face_dir_0to1,
                                       int bc_model_f,
                                       double bc_value_f,
                                       double* pres,
                                       double* dkdp_cell,
                                       WhetStone::DenseMatrix& Jpp)
{
  const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
  double dKrel_dp[2];
  double dpres;

  if (mcells == 2) {
    dpres = pres[0] - pres[1]; // + grn;
    if (little_k_ == OPERATOR_LITTLE_K_UPWIND) {
      double flux0to1;
      flux0to1 = trans_face[0][f] * dpres;
      if (flux0to1 > OPERATOR_UPWIND_RELATIVE_TOLERANCE) { // Upwind
        dKrel_dp[0] = dkdp_cell[0];
        dKrel_dp[1] = 0.0;
      } else if (flux0to1 < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) { // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dkdp_cell[1];
      } else if (fabs(flux0to1) < OPERATOR_UPWIND_RELATIVE_TOLERANCE) { // Upwind
        dKrel_dp[0] = 0.5 * dkdp_cell[0];
        dKrel_dp[1] = 0.5 * dkdp_cell[1];
      }
    } else if (little_k_ == OPERATOR_UPWIND_ARITHMETIC_AVERAGE) {
      dKrel_dp[0] = 0.5 * dkdp_cell[0];
      dKrel_dp[1] = 0.5 * dkdp_cell[1];
    } else {
      AMANZI_ASSERT(0);
    }

    Jpp(0, 0) = trans_face[0][f] * dpres * dKrel_dp[0];
    Jpp(0, 1) = trans_face[0][f] * dpres * dKrel_dp[1];

    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_model_f == OPERATOR_BC_DIRICHLET) {
      pres[1] = bc_value_f;
      dpres = pres[0] - pres[1];
      Jpp(0, 0) = trans_face[0][f] * dpres * dkdp_cell[0];
    } else {
      Jpp(0, 0) = 0.0;
    }
  }
}


/* ******************************************************************
* Compute transmissibilities on faces
****************************************************************** */
void
PDE_DiffusionFV::ComputeTransmissibility_()
{
  Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);

  // Compute auxiliary structure. Note that first components of both
  // fields are used symmetrically (no specific order).
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  CompositeVector beta(cvs, true);
  Epetra_MultiVector& beta_face = *beta.ViewComponent("face", true);
  beta.PutScalar(0.0);

  AmanziMesh::Entity_ID_View faces;
  AmanziGeometry::Point a_dist;
  WhetStone::Tensor Kc(mesh_->getSpaceDimension(), 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c < ncells_owned; ++c) {
    if (K_.get()) Kc = (*K_)[c];
    auto [faces,bisectors] = mesh_->getCellFacesAndBisectors(c);
    int nfaces = faces.size();

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      const AmanziGeometry::Point& a = bisectors[i];
      const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
      double area = mesh_->getFaceArea(f);

      double h_tmp = norm(a);
      double s = area / h_tmp;
      double perm = ((Kc * a) * normal) * s;
      double dxn = a * normal;
      beta_face[0][f] += fabs(dxn / perm);
    }
  }

  beta.GatherGhostedToMaster(Add);

  // Compute transmissibilities. Since it is done only once, we repeat
  // some calculatons.
  transmissibility_->PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++) { trans_face[0][f] = 1.0 / beta_face[0][f]; }
  transmissibility_->ScatterMasterToGhosted("face", true);
  transmissibility_initialized_ = true;
}


/* ******************************************************************
* Return transmissibility value on the given face f.
****************************************************************** */
double
PDE_DiffusionFV::ComputeTransmissibility(int f) const
{
  const Epetra_MultiVector& trans_face = *transmissibility_->ViewComponent("face", true);
  return trans_face[0][f];
}

} // namespace Operators
} // namespace Amanzi
