/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Base class for testing coupled nonlinear diffusion problems:

    - div( k00(u,v) K00 grad(u) ) - div( k01(u,v) K01 grad(v) ) = f0
    - div( k10(u,v) K10 grad(u) ) - div( k11(u,v) K11 grad(v) ) = f1
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_NONLINEAR_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_NONLINEAR_BASE_HH_

#include "Mesh.hh"

class AnalyticNonlinearCoupledBase {
 public:
  AnalyticNonlinearCoupledBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh){};
  virtual ~AnalyticNonlinearCoupledBase() = default;

  // analytic solution for diffusion problem with gravity
  virtual bool isBlock(int i, int j) = 0;

  // -- diffusion tensor T
  virtual Amanzi::WhetStone::Tensor Tensor00(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }
  virtual Amanzi::WhetStone::Tensor Tensor01(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }
  virtual Amanzi::WhetStone::Tensor Tensor10(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }
  virtual Amanzi::WhetStone::Tensor Tensor11(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K(2, 1);
    K(0, 0) = 1.0;
    return K;
  }

  // -- scalar component of the coefficient
  virtual double ScalarCoefficient00(double u, double v) { return 1.; }
  virtual double ScalarCoefficient01(double u, double v) { return 1.; }
  virtual double ScalarCoefficient10(double u, double v) { return 1.; }
  virtual double ScalarCoefficient11(double u, double v) { return 1.; }

  // -- deriviative of the coefficient
  virtual double DScalarCoefficient00D0(double u, double v) { return 0.; }
  virtual double DScalarCoefficient00D1(double u, double v) { return 0.; }
  virtual double DScalarCoefficient01D0(double u, double v) { return 0.; }
  virtual double DScalarCoefficient01D1(double u, double v) { return 0.; }
  virtual double DScalarCoefficient10D0(double u, double v) { return 0.; }
  virtual double DScalarCoefficient10D1(double u, double v) { return 0.; }
  virtual double DScalarCoefficient11D0(double u, double v) { return 0.; }
  virtual double DScalarCoefficient11D1(double u, double v) { return 0.; }

  // -- analytic solution p
  virtual double exact0(const Amanzi::AmanziGeometry::Point& p, double t) const = 0;
  virtual double exact1(const Amanzi::AmanziGeometry::Point& p, double t) const = 0;

  // -- gradient of continuous velocity grad(h), where h = p + g z
  virtual Amanzi::AmanziGeometry::Point
  gradient_exact0(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual Amanzi::AmanziGeometry::Point
  gradient_exact1(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- source term
  virtual double source_exact0(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual double source_exact1(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // derived quantity: Darcy velocity -K * grad(h)
  virtual Amanzi::AmanziGeometry::Point
  velocity_exact0(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K = Tensor00(p, t);
    Amanzi::AmanziGeometry::Point g = gradient_exact0(p, t);
    double u = exact0(p, t);
    double v = exact1(p, t);
    double kr = ScalarCoefficient00(u, v);
    return -(K * g) * kr;
  }
  virtual Amanzi::AmanziGeometry::Point
  velocity_exact1(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::WhetStone::Tensor K = Tensor11(p, t);
    Amanzi::AmanziGeometry::Point g = gradient_exact1(p, t);
    double u = exact0(p, t);
    double v = exact1(p, t);
    double kr = ScalarCoefficient11(u, v);
    Amanzi::AmanziGeometry::Point q = -(K * g) * kr;
    return q;
  }

  // error calculation
  void ComputeCellError(Epetra_MultiVector& u,
                        Epetra_MultiVector& v,
                        double t,
                        double& pnorm,
                        double& l2_err,
                        double& inf_err)
  {
    pnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int ncells =
      mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c < ncells; c++) {
      const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      double u_tmp = exact0(xc, t);
      double v_tmp = exact1(xc, t);
      double volume = mesh_->cell_volume(c);

      // std::cout << c << " " << tmp << " " << p[0][c] << std::endl;
      l2_err += std::pow(u_tmp - u[0][c], 2.0) * volume;
      l2_err += std::pow(v_tmp - v[0][c], 2.0) * volume;
      inf_err = std::max(inf_err, fabs(u_tmp - u[0][c]));
      inf_err = std::max(inf_err, fabs(v_tmp - v[0][c]));
      pnorm += std::pow(u_tmp, 2.0) * volume;
      pnorm += std::pow(v_tmp, 2.0) * volume;
    }
#ifdef HAVE_MPI
    double tmp = pnorm;
    mesh_->get_comm()->SumAll(&tmp, &pnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    pnorm = sqrt(pnorm);
    l2_err = sqrt(l2_err);
  }

  void ComputeFaceError(Epetra_MultiVector& qu,
                        Epetra_MultiVector& qv,
                        double t,
                        double& qnorm,
                        double& l2_err,
                        double& inf_err)
  {
    qnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int nfaces =
      mesh_->num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED);
    for (int f = 0; f < nfaces; f++) {
      double area = mesh_->face_area(f);
      const Amanzi::AmanziGeometry::Point& normal = mesh_->face_normal(f);
      const Amanzi::AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const Amanzi::AmanziGeometry::Point& u_velocity = velocity_exact0(xf, t);
      const Amanzi::AmanziGeometry::Point& v_velocity = velocity_exact1(xf, t);
      double u_tmp = u_velocity * normal;
      double v_tmp = v_velocity * normal;

      l2_err += std::pow((u_tmp - qu[0][f]) / area, 2.0);
      l2_err += std::pow((v_tmp - qv[0][f]) / area, 2.0);
      inf_err = std::max(inf_err, fabs(u_tmp - qu[0][f]) / area);
      inf_err = std::max(inf_err, fabs(v_tmp - qv[0][f]) / area);
      qnorm += std::pow(u_tmp / area, 2.0);
      qnorm += std::pow(v_tmp / area, 2.0);
      // std::cout << f << " " << u[0][f] << " " << tmp << std::endl;
    }
#ifdef HAVE_MPI
    double tmp = qnorm;
    mesh_->get_comm()->SumAll(&tmp, &qnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    qnorm = sqrt(qnorm);
    l2_err = sqrt(l2_err);
  }

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
};

#endif
