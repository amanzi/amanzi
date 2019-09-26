/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for remap methods.
*/

#ifndef AMANZI_OPERATOR_REMAP_DG_TESTS_HH_
#define AMANZI_OPERATOR_REMAP_DG_TESTS_HH_

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "RemapDG.hh"

#include "MeshDeformation.hh"

namespace Amanzi {

template<class AnalyticDG>
class RemapDG_Tests : public Operators::RemapDG {
 public:
  RemapDG_Tests(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                Teuchos::ParameterList& plist) 
    : RemapDG(mesh0, mesh1, plist),
      tprint_(0.0),
      l2norm_(-1.0),
      dt_output_(0.1) {};
  ~RemapDG_Tests() {};

  // experimental options
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
template<class AnalyticDG>
void RemapDG_Tests<AnalyticDG>::InitializeConsistentJacobianDeterminant()
{
  int det_method_tmp = det_method_;
  det_method_ = Operators::OPERATOR_DETERMINANT_VEM;

  // constant part of determinant
  op_adv_->Setup(velc_, false);
  op_adv_->UpdateMatrices(0.0);

  op_reac_->Setup(det_, false);
  op_reac_->UpdateMatrices(0.0);

  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

  op_flux_->Setup(velf_.ptr(), false);
  op_flux_->UpdateMatrices(0.0);
  op_flux_->ApplyBCs(true, true, true);

  CompositeVector& tmp = *op_reac_->global_operator()->rhs();
  CompositeVector one(tmp), u0(tmp), u1(tmp);
  Epetra_MultiVector& one_c = *one.ViewComponent("cell", true);

  one.PutScalarMasterAndGhosted(0.0);
  for (int c = 0; c < ncells_wghost_; ++c) one_c[0][c] = 1.0;

  op_flux_->global_operator()->Apply(one, tmp);
  op_reac_->global_operator()->Apply(tmp, u0);

  // linear part of determinant
  double dt(0.01);
  op_adv_->UpdateMatrices(dt);

  op_flux_->UpdateMatrices(dt);
  op_flux_->ApplyBCs(true, true, true);

  op_flux_->global_operator()->Apply(one, tmp);
  op_reac_->global_operator()->Apply(tmp, u1);
  u1.Update(-1.0/dt, u0, 1.0/dt);

  // save as polynomials
  int nk = one_c.NumVectors();
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
* Print statistics using conservative field u
***************************************************************** */
template<class AnalyticDG>
void RemapDG_Tests<AnalyticDG>::CollectStatistics(double t, const CompositeVector& u)
{
  double tglob = global_time(t);
  if (tglob >= tprint_) {
    op_reac_->UpdateMatrices(Teuchos::null);
    auto& matrices = op_reac_->local_op()->matrices;
    for (int n = 0; n < matrices.size(); ++n) matrices[n].Inverse();

    auto& rhs = *op_reac_->global_operator()->rhs();
    op_reac_->global_operator()->Apply(u, rhs);
    rhs.Dot(u, &l2norm_);

    Epetra_MultiVector& xc = *rhs.ViewComponent("cell");
    int nk = xc.NumVectors();
    double xmax[nk], xmin[nk];
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (mesh0_->get_comm()->MyPID() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  umax: ", tglob, l2norm_, nfun_, sharp_);
      for (int i = 0; i < std::min(nk, 4); ++i) printf("%9.5g ", xmax[i]);
      printf("  umin: %9.5g\n", xmin[0]);
    }

    tprint_ += dt_output_;
    sharp_ = 0.0;
  } 
}

} // namespace Amanzi

#endif
