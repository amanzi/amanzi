/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Base class for remap methods.
*/

#ifndef AMANZI_OPERATOR_REMAP_DG_TESTS_HH_
#define AMANZI_OPERATOR_REMAP_DG_TESTS_HH_

// Amanzi::Operators
#include "OperatorDefs.hh"
#include "RemapDG.hh"

#include "MeshDeformation.hh"

namespace Amanzi {

template <class AnalyticDG>
class MyRemapDGBase : public Operators::RemapDG {
 public:
  MyRemapDGBase(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                Teuchos::ParameterList& plist)
    : RemapDG(mesh0, mesh1, plist), tprint_(0.0), l2norm_(-1.0), dt_output_(0.1){};
  ~MyRemapDGBase(){};

  // CFL condition
  double StabilityCondition();

  // tools
  // -- mass on mesh0
  double InitialMass(const CompositeVector& p1, int order);

  // output
  void CollectStatistics(double t, const CompositeVector& u);
  virtual double global_time(double t) { return t; }
  void set_dt_output(double dt) { dt_output_ = dt; }

  // statistics
  double tprint_, dt_output_, l2norm_;
};


/* *****************************************************************
* Initialization of the consistent jacobian determinant
***************************************************************** */
template <class AnalyticDG>
double
MyRemapDGBase<AnalyticDG>::StabilityCondition()
{
  double dt(1e+99), alpha(0.2), tmp;

  for (int f = 0; f < nfaces_wghost_; ++f) {
    double area = mesh0_->face_area(f);
    const AmanziGeometry::Point& xf = mesh0_->face_centroid(f);
    velf_vec_[f].Value(xf).Norm2(&tmp);
    dt = std::min(dt, area / tmp);
  }

  return dt * alpha / (2 * order_ + 1);
}


/* *****************************************************************
* Compute initial mass
***************************************************************** */
template <class AnalyticDG>
double
MyRemapDGBase<AnalyticDG>::InitialMass(const CompositeVector& p1, int order)
{
  const Epetra_MultiVector& p1c = *p1.ViewComponent("cell", false);
  int nk = p1c.NumVectors();
  int ncells = p1c.MyLength();

  double mass(0.0), mass0;
  WhetStone::DenseVector data(nk);
  WhetStone::NumericalIntegration<AmanziMesh::Mesh> numi(mesh0_);

  for (int c = 0; c < ncells; c++) {
    for (int i = 0; i < nk; ++i) data(i) = p1c[i][c];
    auto poly = dg_->cell_basis(c).CalculatePolynomial(mesh0_, c, order, data);
    mass += numi.IntegratePolynomialCell(c, poly);
  }

  mesh0_->get_comm()->SumAll(&mass, &mass0, 1);
  return mass0;
}


/* *****************************************************************
* Print statistics using conservative field u
***************************************************************** */
template <class AnalyticDG>
void
MyRemapDGBase<AnalyticDG>::CollectStatistics(double t, const CompositeVector& u)
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
    double xmax[nk], xmin[nk], lmax(-1.0), lmin(-1.0), lavg(-1.0);
    xc.MaxValue(xmax);
    xc.MinValue(xmin);

    if (limiter() != Teuchos::null) {
      const auto& lim = *limiter()->limiter();
      lim.MaxValue(&lmax);
      lim.MinValue(&lmin);
      lim.MeanValue(&lavg);
    }

    if (mesh0_->get_comm()->MyPID() == 0) {
      printf("t=%8.5f  L2=%9.5g  nfnc=%5d  sharp=%5.1f%%  limiter: %6.3f %6.3f %6.3f  umax/umin: "
             "%9.5g %9.5g\n",
             tglob,
             l2norm_,
             nfun_,
             sharp_,
             lmax,
             lmin,
             lavg,
             xmax[0],
             xmin[0]);
    }

    tprint_ += dt_output_;
    sharp_ = 0.0;
  }
}

} // namespace Amanzi

#endif
