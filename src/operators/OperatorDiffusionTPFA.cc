/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "Epetra_FECrsGraph.h"

#include "LinearOperatorFactory.hh"
#include "mfd3d_diffusion.hh"
#include "PreconditionerFactory.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionTPFA.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Constructor.                                           
****************************************************************** */
void OperatorDiffusionTPFA::InitOperator(
    std::vector<WhetStone::Tensor>& K,
    Teuchos::RCP<const CompositeVector> k, Teuchos::RCP<const CompositeVector> dkdp,
    double rho, double mu)
{
  K_ = &K;
  k_ = k;
  dkdp_ = dkdp;

  rho_ = rho;
  mu_ = mu;
  scalar_rho_mu_ = true;

  // solution-independent data
  const Epetra_BlockMap& fmap_wghost = mesh_->face_map(true);
  transmissibility_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  gravity_term_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));

  ComputeTransmissibilities_();
}


/* ******************************************************************
* Populate face-based matrices.
****************************************************************** */
void OperatorDiffusionTPFA::UpdateMatrices(Teuchos::RCP<const CompositeVector> flux, 
                                           Teuchos::RCP<const CompositeVector> p)
{
  const std::vector<int>& bc_model = bc_->bc_model();

  // find location of matrix blocks
  int schema_my = OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m(0), nblocks = blocks_.size();
  bool flag(false);

  for (int n = 0; n < nblocks; n++) {
    if (blocks_schema_[n] == schema_my) {
      m = n;
      flag = true;
      break;
    }
  }

  if (flag == false) { 
    m = nblocks++;
    blocks_schema_.push_back(OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_DOFS_CELL);
    blocks_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
    blocks_shadow_.push_back(Teuchos::rcp(new std::vector<WhetStone::DenseMatrix>));
  }
  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];
  WhetStone::DenseMatrix null_matrix;

  // preparing upwind data
  const Epetra_MultiVector& k_face = *k_->ViewComponent("face", true);

  // updating matrix blocks
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    WhetStone::DenseMatrix Aface(ncells, ncells);

    double tij = (*transmissibility_)[f] * k_face[0][f];
    for (int i = 0; i < ncells; i++) {
      Aface(i, i) = tij;
      for (int j = i + 1; j < ncells; j++) {
	Aface(i, j) = -tij;
	Aface(j, i) = -tij;
      }
    }

    if (flag) {
      matrix[f] += Aface;
    } else {
      matrix.push_back(Aface);
      matrix_shadow.push_back(null_matrix);
    }
  }

  // populating right-hand side
  Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      if (bc_model[f] == OPERATOR_BC_FACE_NEUMANN) continue;
      rhs_cell[0][c] -= dirs[n] * (*gravity_term_)[f] * k_face[0][f];  
    }
  }

  // Add derivatives to the matrix (Jacobian in this case)
  AnalyticJacobian_(*p);
}


/* ******************************************************************
* Special implementation of boundary conditions.
****************************************************************** */
void OperatorDiffusionTPFA::ApplyBCs()
{
  const std::vector<int>& bc_model = bc_->bc_model();
  const std::vector<double>& bc_value = bc_->bc_value();

  const Epetra_MultiVector& k_face = *k_->ViewComponent("face", true);
  Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");

  AmanziMesh::Entity_ID_List cells;

  for (int f = 0; f < nfaces_wghost; f++) {
    if (bc_model[f] != OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];

      if (bc_model[f] == OPERATOR_BC_FACE_DIRICHLET) {
        rhs_cell[0][c] += bc_value[f] * (*transmissibility_)[f] * k_face[0][f];
      } else if (bc_model[f] == OPERATOR_BC_FACE_NEUMANN) {
        rhs_cell[0][c] -= bc_value[f] * mesh_->face_area(f);
        (*transmissibility_)[f] = 0.0;
        (*gravity_term_)[f] = 0.0;
      }
    }
  }
}


/* ******************************************************************
* Parallel matvec product Y = A_ * X.
****************************************************************** */
int OperatorDiffusionTPFA::Apply(const CompositeVector& X, CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

  int ierr = A_->Multiply(false, Xc, Yc);

  if (ierr) {
    Errors::Message msg("OperatorDiffusionTPFA::Apply has failed to calculate y = A*x.");
    Exceptions::amanzi_throw(msg);
  }

  return 0;
}


/* ******************************************************************
* Invert Jacobian with a GMRES solver.
****************************************************************** */
int OperatorDiffusionTPFA::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& pre_list = plist.sublist("gmres");
  Teuchos::ParameterList& slist = pre_list.sublist("gmres parameters");

  pre_list.set<string>("iterative method", "gmres");
  slist.set<double>("error tolerance", 1e-7);
  slist.set<int>("maximum number of iterations", 200);
  Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
  vlist.set("Verbosity Level", "low");

  // delegating preconditioning to the base operator
  Teuchos::RCP<const OperatorDiffusionTPFA> op_matrix = Teuchos::rcp(this, false);
  Teuchos::RCP<const Operator> op_preconditioner = Teuchos::rcp(new Operator(*this));

  AmanziSolvers::LinearOperatorFactory<Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operator, CompositeVector, CompositeVectorSpace> > 
      solver = factory.Create("gmres", plist, op_matrix, op_preconditioner);

  Y.PutScalar(0.0);
  int ierr = solver->ApplyInverse(X, Y);

  return 0;
}


/* ******************************************************************
* Linear algebra operations with matrices: r = A * u - f                                                 
****************************************************************** */
void OperatorDiffusionTPFA::ComputeNegativeResidual(const CompositeVector& u, CompositeVector& r)
{
  const Epetra_MultiVector& k_face = *k_->ViewComponent("face", true);
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);
  Epetra_MultiVector& rc = *r.ViewComponent("cell");

  AmanziMesh::Entity_ID_List cells;

  // matvec product A*u
  u.ScatterMasterToGhosted("cell");

  r.PutScalar(0.0);

  for (int f = 0; f < nfaces_wghost; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      int c1 = cells[0];
      int c2 = cells[1];

      double tmp = k_face[0][f] * (*transmissibility_)[f] * (uc[0][c1] - uc[0][c2]);  
      if (c1 < ncells_owned) rc[0][c1] += tmp; 
      if (c2 < ncells_owned) rc[0][c2] -= tmp; 
    } else {
      int c = cells[0];
      double tmp = k_face[0][f] * (*transmissibility_)[f] * uc[0][c];
      if (c < ncells_owned) rc[0][c] += tmp;
    }							
  } 
  
  // right-hand side
  Epetra_MultiVector& rhs_cell = *rhs_->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {    
    rc[0][c] -= rhs_cell[0][c];   
  }
}


/* ******************************************************************
* Calculate flux from cell-centered data
****************************************************************** */
void OperatorDiffusionTPFA::UpdateFlux(
    const CompositeVector& solution, CompositeVector& darcy_mass_flux)
{
  const std::vector<int>& bc_model = bc_->bc_model();
  const std::vector<double>& bc_value = bc_->bc_value();

  solution.ScatterMasterToGhosted("cell");

  const Epetra_MultiVector& Krel_face = *k_->ViewComponent("face", true);
  const Epetra_MultiVector& p = *solution.ViewComponent("cell", true);
  Epetra_MultiVector& flux = *darcy_mass_flux.ViewComponent("face", true);

  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  std::vector<int> flag(nfaces_wghost, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] == OPERATOR_BC_FACE_DIRICHLET) {
	double value = bc_value[f];
	flux[0][f] = dirs[n] * (*transmissibility_)[f] * (p[0][c] - value) + (*gravity_term_)[f];
	flux[0][f] *= Krel_face[0][f];

      } else if (bc_model[f] == OPERATOR_BC_FACE_NEUMANN) {
	double value = bc_value[f];
	double area = mesh_->face_area(f);
	flux[0][f] = value*area;

      } else {
	if (f < nfaces_owned && !flag[f]) {
	  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
	  if (cells.size() <= 1) {
	    Errors::Message msg("Flow PK: These boundary conditions are not supported by TPFA.");
	    Exceptions::amanzi_throw(msg);
	  }
          int c1 = cells[0];
          int c2 = cells[1];
          if (c == c1) {
            flux[0][f] = dirs[n] * (*transmissibility_)[f] * (p[0][c1] - p[0][c2]) + (*gravity_term_)[f];
          } else {
            flux[0][f] = dirs[n] * (*transmissibility_)[f] * (p[0][c2] - p[0][c1]) + (*gravity_term_)[f];
          }	    
          flux[0][f] *= Krel_face[0][f];
          flag[f] = 1;
	}
      }
    }
  }
}


/* ******************************************************************
* Computation the part of the Jacobian which depends on derivatives 
* of the relative permeability wrt to capillary pressure. They must
* be added to the existing matrix structure.
****************************************************************** */
void OperatorDiffusionTPFA::AnalyticJacobian_(const CompositeVector& u)
{
  const std::vector<int>& bc_model = bc_->bc_model();
  const std::vector<double>& bc_value = bc_->bc_value();

  // find location of matrix blocks
  int schema_my = OPERATOR_SCHEMA_BASE_FACE + OPERATOR_SCHEMA_DOFS_CELL;
  int m(0), nblocks = blocks_.size();

  for (int n = 0; n < nblocks; n++) {
    if (blocks_schema_[n] == schema_my) {
      m = n;
      break;
    }
  }

  std::vector<WhetStone::DenseMatrix>& matrix = *blocks_[m];
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = *blocks_shadow_[m];

  u.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& uc = *u.ViewComponent("cell", true);

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  AmanziMesh::Entity_ID_List cells, faces;
  std::vector<int> dirs;

  int gid[2];
  double k_rel[2], dkdp[2], pres[2], dist;

  const Epetra_MultiVector& dKdP_cell = *dkdp_->ViewComponent("cell");
  const Epetra_MultiVector& dKdP_face = *dkdp_->ViewComponent("face", true);

  std::vector<int> flag(nfaces_owned, 0);

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int i = 0; i < nfaces; i++) {
      int f = faces[i];
      if (f < nfaces_owned && !flag[f]) {
	mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
	int mcells = cells.size();

	WhetStone::DenseMatrix Aface(mcells, mcells);

	int face_dir;   	
	// face_dir is equal to 1 if normal direction points from cells[0],
	// otherwise face_dir is -1
	if (cells[0] == c) face_dir = dirs[i];
	else face_dir = -dirs[i];

	for (int n = 0; n < mcells; n++) {
          int c1 = cells[n];
	  gid[n] = cmap_wghost.GID(c1);
	  pres[n] = uc[0][c1];
	  dkdp[n] = dKdP_cell[0][c1];
	}

	if (mcells == 1) dkdp[0] = dKdP_face[0][f];

	const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, cells[0]);
	ComputeJacobianLocal_(mcells, f, face_dir, upwind_, bc_model[f], bc_value[f], pres, dkdp, Aface);

        matrix[f] += Aface;
	flag[f] = 1;
      }
    }
  }
}


/* ******************************************************************
* Computation of a local submatrix of the analytical Jacobian 
* (its nonlinear part) on face f.
****************************************************************** */
void OperatorDiffusionTPFA::ComputeJacobianLocal_(
    int mcells, int f, int face_dir, int Krel_method,
    int bc_models, double bc_value,
    double *pres, double *dkdp_cell,
    WhetStone::DenseMatrix& Jpp)
{
  double dKrel_dp[2];

  double dpres;
  if (mcells == 2) {
    dpres = pres[0] - pres[1];  // + grn;
    if (Krel_method == OPERATOR_UPWIND_WITH_CONSTANT_VECTOR) {  // Define K 
      double cos_angle = face_dir * (*gravity_term_)[f] / (rho_ * mesh_->face_area(f));
      if (cos_angle > OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dkdp_cell[0];
        dKrel_dp[1] = 0.0;
      } else if (cos_angle < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dkdp_cell[1];
      } else if (fabs(cos_angle) < OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.5 * dkdp_cell[0];
        dKrel_dp[1] = 0.5 * dkdp_cell[1];
      }
    } else if (Krel_method == OPERATOR_UPWIND_WITH_FLUX) {
      double flux0to1;
      flux0to1 = (*transmissibility_)[f] * dpres + face_dir * (*gravity_term_)[f];
      if (flux0to1  > OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = dkdp_cell[0];
        dKrel_dp[1] = 0.0;
      } else if ( flux0to1 < -OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.0;
        dKrel_dp[1] = dkdp_cell[1];
      } else if (fabs(flux0to1) < OPERATOR_UPWIND_RELATIVE_TOLERANCE) {  // Upwind
        dKrel_dp[0] = 0.5 * dkdp_cell[0];
        dKrel_dp[1] = 0.5 * dkdp_cell[1];
      }
    } else if (Krel_method == OPERATOR_ARITHMETIC_MEAN) {
      dKrel_dp[0] = 0.5 * dkdp_cell[0];
      dKrel_dp[1] = 0.5 * dkdp_cell[1];
    }

    Jpp(0, 0) = ((*transmissibility_)[f] * dpres + face_dir * (*gravity_term_)[f]) * dKrel_dp[0];
    Jpp(0, 1) = ((*transmissibility_)[f] * dpres + face_dir * (*gravity_term_)[f]) * dKrel_dp[1];
    Jpp(1, 0) = -Jpp(0, 0);
    Jpp(1, 1) = -Jpp(0, 1);

  } else if (mcells == 1) {
    if (bc_models == OPERATOR_BC_FACE_DIRICHLET) {                   
      pres[1] = bc_value;
      dpres = pres[0] - pres[1];  // + grn;
      Jpp(0, 0) = ((*transmissibility_)[f] * dpres + face_dir * (*gravity_term_)[f]) * dkdp_cell[0];
    } else {
      Jpp(0, 0) = 0.0;
    }
  }
}


/* ******************************************************************
* Compute transmissibilities on faces 
****************************************************************** */
void OperatorDiffusionTPFA::ComputeTransmissibilities_()
{
  transmissibility_->PutScalar(0.0);

  AmanziGeometry::Point gravity(g_ * rho_);

  AmanziMesh::Entity_ID_List cells;
  AmanziGeometry::Point a_dist, a[2];
  double h[2], perm[2], beta[2], trans_f;

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    const AmanziGeometry::Point& normal = mesh_->face_normal(f, false);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    if (ncells == 2) {
      a_dist = mesh_->cell_centroid(cells[1]) - mesh_->cell_centroid(cells[0]);
    } else if (ncells == 1) {    
      a_dist = xf - mesh_->cell_centroid(cells[0]);
    } 

    a_dist *= 1.0 / norm(a_dist);

    for (int i = 0; i < ncells; i++) {
      int c = cells[i];
      a[i] = xf - mesh_->cell_centroid(c);
      h[i] = norm(a[i]);
      double s = area / h[i];
      perm[i] = (rho_ / mu_) * (((*K_)[c] * a[i]) * normal) * s;

      double dxn = a[i] * normal;
      beta[i] = fabs(perm[i] / dxn);
    }

    double grav = (gravity * normal) / area;
    trans_f = 0.0;

    if (ncells == 2) {
      grav *= (h[0] + h[1]);
      trans_f = (beta[0] * beta[1]) / (beta[0] + beta[1]);
    } else if (ncells == 1) {    
      grav *= h[0];
      trans_f = beta[0];
    } 

    (*transmissibility_)[f] = trans_f;
    (*gravity_term_)[f] = (*transmissibility_)[f] * grav;
  }

  // parallelization using CV capability
#ifdef HAVE_MPI
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("face", AmanziMesh::FACE, 1);

  CompositeVector tmp(cvs, true);
  Epetra_MultiVector& data = *tmp.ViewComponent("face", true);

  data = *transmissibility_;
  tmp.ScatterMasterToGhosted("face", true);
  for (int f = nfaces_owned; f < nfaces_wghost; f++) {
    (*transmissibility_)[f] = data[0][f];
  }

  data = *gravity_term_;
  tmp.ScatterMasterToGhosted("face", true);
  for (int f = nfaces_owned; f < nfaces_wghost; f++) {
    (*gravity_term_)[f] = data[0][f];
  }
#endif
}

}  // namespace Operators
}  // namespace Amanzi


