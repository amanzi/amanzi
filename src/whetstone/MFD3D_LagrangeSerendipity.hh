/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Serendipity Lagrange-type element: degrees of freedom are nodal values 
  and moments on edges, faces and inside cell. The number of later is 
  reduced significantly for polygonal cells. 
*/

#ifndef AMANZI_MFD3D_LAGRANGE_SERENDIPITY_HH_
#define AMANZI_MFD3D_LAGRANGE_SERENDIPITY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MFD3D_Lagrange.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeSerendipity : public MFD3D_Lagrange { 
 public:
  MFD3D_LagrangeSerendipity(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : InnerProduct(mesh),
      MFD3D_Lagrange(mesh) {};
  ~MFD3D_LagrangeSerendipity() {};

  // required methods
  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- projectors
  virtual void L2Cell(
      int c, const std::vector<Polynomial>& vf,
      Polynomial& moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, Type::L2, moments, uc);
  }

  virtual void H1Cell(
      int c, const std::vector<Polynomial>& vf,
      Polynomial& moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, Type::H1, moments, uc);
  }

  // other methods
  void L2Cell_LeastSquare(
      int c, const std::vector<Polynomial>& vf,
      Polynomial& moments, Polynomial& uc) {
    ProjectorCell_(c, vf, Type::LS, moments, uc);
  }

 private:
  void ProjectorCell_(
      int c, const std::vector<Polynomial>& vf,
      const Projectors::Type type,
      Polynomial& moments, Polynomial& uc);

  void CalculateDOFsOnBoundary_(
      int c, const std::vector<Polynomial>& vf, DenseVector& vdof);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

