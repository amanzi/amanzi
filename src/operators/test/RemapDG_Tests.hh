/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_REMAP_DG_TESTS_HH_
#define AMANZI_OPERATOR_REMAP_DG_TESTS_HH_

#include "OperatorDefs.hh"
#include "RemapDG.hh"

namespace Amanzi {

template <class AnalyticDG>
class RemapDG_Tests : public Operators::RemapDG {
 public:
  RemapDG_Tests(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                Teuchos::ParameterList& plist)
    : RemapDG(mesh0, mesh1, plist),
      tprint_(0.0),
      l2norm_(-1.0),
      dt_output_(0.1){};
  ~RemapDG_Tests(){};

  // mesh deformation
  virtual void DeformMesh(int deform, double t);
  AmanziGeometry::Point
  DeformNode(int deform, double t, const AmanziGeometry::Point& yv,
             const AmanziGeometry::Point& rv = AmanziGeometry::Point(3));

  // experimental options
  virtual void DynamicCellVelocity(double t);
  void InitializeConsistentJacobianDeterminant();

  // output
  void CollectStatistics(double t, const CompositeVector& u);
  virtual double global_time(double t) { return t; }
  void set_dt_output(double dt) { dt_output_ = dt; }

 protected:
  std::vector<WhetStone::Polynomial> det0_, det1_;

  // statistics
  double tprint_, dt_output_, l2norm_;
};


/* *****************************************************************
 * Initialization of the consistent jacobian determinant
 ***************************************************************** */
template <class AnalyticDG>
void
RemapDG_Tests<AnalyticDG>::InitializeConsistentJacobianDeterminant()
{
  int det_method_tmp = det_method_;
  det_method_ = Operators::OPERATOR_DETERMINANT_VEM;

  // constant part of determinant
  DynamicFaceVelocity(0.0);
  DynamicCellVelocity(0.0);

  op_adv_->SetupPolyVector(velc_);
  op_adv_->UpdateMatrices();

  op_reac_->Setup(det_);
  op_reac_->UpdateMatrices(Teuchos::null);

  auto& matrices = op_reac_->local_matrices()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

  op_flux_->Setup(velc_, velf_);
  op_flux_->UpdateMatrices(velf_.ptr());
  op_flux_->ApplyBCs(true, true, true);

  CompositeVector& tmp = *op_reac_->global_operator()->rhs();
  CompositeVector one(tmp), u0(tmp), u1(tmp);
  Epetra_MultiVector& one_c = *one.ViewComponent("cell", true);

  one.putScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_wghost_; ++c) one_c[0][c] = 1.0;

  op_flux_->global_operator()->apply(one, tmp);
  op_reac_->global_operator()->apply(tmp, u0);

  // linear part of determinant
  double dt(0.01);
  DynamicFaceVelocity(dt);
  DynamicCellVelocity(dt);

  op_adv_->SetupPolyVector(velc_);
  op_adv_->UpdateMatrices();

  op_flux_->Setup(velc_, velf_);
  op_flux_->UpdateMatrices(velf_.ptr());
  op_flux_->ApplyBCs(true, true, true);

  op_flux_->global_operator()->apply(one, tmp);
  op_reac_->global_operator()->apply(tmp, u1);
  u1.update(-1.0 / dt, u0, 1.0 / dt);

  // save as polynomials
  int nk = one_c.getNumVectors();
  Amanzi::WhetStone::DenseVector data(nk);
  Epetra_MultiVector& u0c = *u0.ViewComponent("cell", true);
  Epetra_MultiVector& u1c = *u1.ViewComponent("cell", true);

  det0_.resize(ncells_owned_);
  det1_.resize(ncells_owned_);

  for (int c = 0; c < ncells_owned_; ++c) {
    const auto& basis = dg_->cell_basis(c);

    for (int i = 0; i < nk; ++i) data(i) = u0c[i][c];
    det0_[c] = basis.CalculatePolynomial(mesh0_, c, order_, data);

    for (int i = 0; i < nk; ++i) data(i) = u1c[i][c];
    det1_[c] = basis.CalculatePolynomial(mesh0_, c, order_, data);
  }

  det_method_ = det_method_tmp;
}


/* *****************************************************************
 * Cell co-velocity in reference coordinates and Jacobian determinant
 ***************************************************************** */
template <class AnalyticDG>
void
RemapDG_Tests<AnalyticDG>::DynamicCellVelocity(double t)
{
  WhetStone::MatrixPolynomial Jt, C;
  for (int c = 0; c < ncells_owned_; ++c) {
    DynamicJacobianMatrix(c, t, J_[c], Jt);
    maps_->Cofactors(Jt, C);
    if (det_method_ == Operators::OPERATOR_DETERMINANT_EXACT_TI) {
      double tmp = t * t / 2;
      (*det_)[c][0] = t * det0_[c] + tmp * det1_[c];
      (*det_)[c][0](0) += 1.0;
    } else if (det_method_ == Operators::OPERATOR_DETERMINANT_VEM) {
      maps_->Determinant(Jt, (*det_)[c]);
    }

    // negative co-velocity, v = -C^t u
    int nC = C.NumRows();
    (*velc_)[c].resize(nC);

    int kC = nC / dim_;
    for (int n = 0; n < kC; ++n) {
      int m = n * dim_;
      for (int i = 0; i < dim_; ++i) {
        (*velc_)[c][m + i].Reshape(dim_, 0, true);
        (*velc_)[c][m + i].set_origin(uc_[c][0].origin());

        for (int k = 0; k < dim_; ++k) {
          (*velc_)[c][m + i] -= C(m + k, i) * uc_[c][m + k];
        }
      }
    }
  }
}


/* *****************************************************************
 * Print statistics using conservative field u
 ***************************************************************** */
template <class AnalyticDG>
void
RemapDG_Tests<AnalyticDG>::CollectStatistics(double t, const CompositeVector& u)
{
  double tglob = global_time(t);
  if (tglob >= tprint_) {
    op_reac_->UpdateMatrices(Teuchos::null);
    auto& matrices = op_reac_->local_matrices()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();

    auto& rhs = *op_reac_->global_operator()->rhs();
    op_reac_->global_operator()->apply(u, rhs);
    l2norm_ = rhs.dot(u);

    Epetra_MultiVector& xc = *rhs.ViewComponent("cell");
    int nk = xc.getNumVectors();
    double xmax[nk], xmin[nk];
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (mesh0_->get_comm()->getRank() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  umax: ",
             tglob,
             l2norm_,
             nfun_,
             sharp_);
      for (int i = 0; i < std::min(nk, 4); ++i) printf("%9.5g ", xmax[i]);
      printf("  umin: %9.5g\n", xmin[0]);
    }

    tprint_ += dt_output_;
    sharp_ = 0.0;
  }
}


/* *****************************************************************
 * Deform mesh1
 ***************************************************************** */
template <class AnalyticDG>
void
RemapDG_Tests<AnalyticDG>::DeformMesh(int deform, double t)
{
  // create distributed random vector
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh0_)->SetGhosted(true)->AddComponent(
    "node", AmanziMesh::NODE, dim_);
  CompositeVector random(cvs);

  int gid = mesh0_->node_map(false).MaxAllGID();
  double scale = 0.2 * std::pow(gid, -1.0 / dim_);
  Epetra_MultiVector& random_n = *random.ViewComponent("node", true);

  random_n.Random();
  random_n.scale(scale);
  random.ScatterMasterToGhosted();

  // relocate mesh nodes
  AmanziGeometry::Point xv(dim_), yv(dim_), uv(dim_), rv(dim_);
  AmanziMesh::Entity_ID_List nodeids;
  AmanziGeometry::Point_List new_positions, final_positions;

  int nnodes =
    mesh0_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  for (int v = 0; v < nnodes; ++v) {
    for (int i = 0; i < dim_; ++i) rv[i] = random_n[i][v];
    mesh0_->node_get_coordinates(v, &xv);
    nodeids.push_back(v);
    new_positions.push_back(DeformNode(deform, t, xv, rv));
  }
  mesh1_->deform(nodeids, new_positions, false, &final_positions);
}


/* *****************************************************************
 * Deformation functional
 ***************************************************************** */
template <class AnalyticDG>
AmanziGeometry::Point
RemapDG_Tests<AnalyticDG>::DeformNode(int deform, double t,
                                      const AmanziGeometry::Point& yv,
                                      const AmanziGeometry::Point& rv)
{
  AmanziGeometry::Point uv(dim_), xv(yv);

  if (deform == 1) {
    double ds(0.0001);
    int n = t / ds;
    for (int i = 0; i < n; ++i) {
      if (dim_ == 2) {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]);
        uv[1] = -0.2 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]);
      } else {
        uv[0] = 0.2 * std::sin(M_PI * xv[0]) * std::cos(M_PI * xv[1]) *
                std::cos(M_PI * xv[2]);
        uv[1] = -0.1 * std::cos(M_PI * xv[0]) * std::sin(M_PI * xv[1]) *
                std::cos(M_PI * xv[2]);
        uv[2] = -0.1 * std::cos(M_PI * xv[0]) * std::cos(M_PI * xv[1]) *
                std::sin(M_PI * xv[2]);
      }
      xv += uv * ds;
    }
  } else if (deform == 2) {
    xv[0] = yv[0] * yv[1] + (1.0 - yv[1]) * std::pow(yv[0], 0.8);
    xv[1] = yv[1] * yv[0] + (1.0 - yv[0]) * std::pow(yv[1], 0.8);
  } else if (deform == 3) {
    if (fabs(yv[0]) > 1e-6 && fabs(1.0 - yv[0]) > 1e-6 && fabs(yv[1]) > 1e-6 &&
        fabs(1.0 - yv[1]) > 1e-6) {
      xv[0] += rv[0];
      xv[1] += rv[1];
    }
  } else if (deform == 4) {
    xv[0] += t * yv[0] * yv[1] * (1.0 - yv[0]) / 2;
    xv[1] += t * yv[0] * yv[1] * (1.0 - yv[1]) / 2;
  } else if (deform == 5) {
    xv[0] += t * yv[0] * (1.0 - yv[0]) / 2;
    xv[1] += t * yv[1] * (1.0 - yv[1]) / 2;
  } else if (deform == 6) {
    double phi = t * 2 * M_PI;
    double cs(std::cos(phi)), sn(std::sin(phi));
    xv[0] = cs * yv[0] - sn * yv[1];
    xv[1] = sn * yv[0] + cs * yv[1];
  }

  return xv;
}

} // namespace Amanzi

#endif
