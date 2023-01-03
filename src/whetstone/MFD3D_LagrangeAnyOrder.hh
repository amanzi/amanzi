/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Lagrange-type element: degrees of freedom are ordered as follows:
    (1) nodal values in the natural order;
    (2) moments on faces groupped by face;
    (3) moments of edges, groupped by edge (in 3D);
    (4) moments inside cell.
*/

#ifndef AMANZI_MFD3D_LAGRANGE_ANY_ORDER_HH_
#define AMANZI_MFD3D_LAGRANGE_ANY_ORDER_HH_

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "GrammMatrix.hh"
#include "MFD3D.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "SurfaceCoordinateSystem.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeAnyOrder : public MFD3D {
 public:
  MFD3D_LagrangeAnyOrder(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : MFD3D(mesh){};
  MFD3D_LagrangeAnyOrder(const Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac)
  {
    if (d_ == 2) return H1consistency2D_(mesh_, c, T, N, Ac);
    return H1consistency3D_(c, T, N, Ac, true);
  }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(c, ve, vf, ProjectorType::L2, moments, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(c, ve, vf, ProjectorType::H1, moments, uc);
  }

  virtual void H1Cell(int c, const DenseVector& dofs, Polynomial& uc) override
  {
    ProjectorCellFromDOFs_(c, dofs, ProjectorType::H1, uc);
  }

  // surface methods
  int StiffnessMatrixSurface(int c, const Tensor& T, DenseMatrix& A);

  // access
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  PolynomialOnMesh& integrals() { return integrals_; }

  // -- matrices that could be resused in other code
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 protected:
  int H1consistency2D_(const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
                       int c,
                       const Tensor& T,
                       DenseMatrix& N,
                       DenseMatrix& Ac);

  int H1consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac, bool doAc);

 private:
  void ProjectorCell_(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments,
                      Polynomial& uc);

  void
  ProjectorCellFromDOFs_(int c, const DenseVector& dofs, const ProjectorType type, Polynomial& uc);

  std::vector<Polynomial> ConvertMomentsToPolynomials_(int order);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_LagrangeAnyOrder> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
