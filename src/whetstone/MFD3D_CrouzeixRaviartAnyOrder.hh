/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#ifndef AMANZI_MFD3D_CROUZEIX_RAVIART_ANY_ORDER_HH_
#define AMANZI_MFD3D_CROUZEIX_RAVIART_ANY_ORDER_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "GrammMatrix.hh"
#include "MatrixObjects.hh"
#include "MFD3D.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"
#include "VectorObjectsUtils.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviartAnyOrder : public MFD3D { 
 public:
  MFD3D_CrouzeixRaviartAnyOrder(const Teuchos::ParameterList& plist,
                                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  void L2GradientCell(int c, const std::vector<VectorPolynomial>& vf,
                      const std::shared_ptr<DenseVector>& moments, MatrixPolynomial& uc) {
    ProjectorGradientCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, ProjectorType::H1, moments, uc);
  }

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  // generic code for multiple projectors
  template<class MyMesh>
  void ProjectorCell_(const Teuchos::RCP<const MyMesh>& mymesh, 
                      int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type, 
                      const Polynomial* moments, Polynomial& uc);

  void ProjectorGradientCell_(int c, const std::vector<VectorPolynomial>& vf,
                              const ProjectorType type, 
                              const std::shared_ptr<DenseVector>& moments, MatrixPolynomial& uc);

  // supporting routines
  void CalculateFaceDOFs_(int f, const Polynomial& vf, const Polynomial& pf,
                          DenseVector& vdof, int& row);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_CrouzeixRaviartAnyOrder> factory_;
};


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
template<class MyMesh>
void MFD3D_CrouzeixRaviartAnyOrder::ProjectorCell_(
    const Teuchos::RCP<const MyMesh>& mymesh, 
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, Polynomial& uc)
{
  AMANZI_ASSERT(d_ == 2);

  Entity_ID_List faces;
  mymesh->cell_get_faces(c, &faces);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mymesh->cell_centroid(c);
  double volume = mymesh->cell_volume(c);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  StiffnessMatrix(c, T, A);  

  // number of degrees of freedom
  Polynomial pf(d_ - 1, order_ - 1);
  int nd = G_.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  // calculate DOFs for boundary polynomial
  DenseVector vdof(ndof);

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized<AmanziMesh::Mesh> basis;
  basis.Init(mymesh, c, order_, ptmp);

  // populate matrices N and R
  int row(0);
  // degrees of freedom on faces
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    CalculateFaceDOFs_(f, vf[n], pf, vdof, row);
  }
  // degrees of freedom in cell
  if (ndof_c > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    AMANZI_ASSERT(ndof_c == v3.NumRows());

    for (int n = 0; n < ndof_c; ++n) {
      vdof(row + n) = v3(n);
    }
  }

  // calculate polynomial coefficients (in natural basis)
  DenseVector v4(nd), v5(nd);
  R_.Multiply(vdof, v4, true);
  G_.Multiply(v4, v5, false);

  uc = basis.CalculatePolynomial(mymesh, c, order_, v5);

  // uniqueness requires to specify constant in polynomial
  if (order_ == 1) {
    AmanziGeometry::Point grad(d_);
    for (int j = 0; j < d_; ++j) {
      grad[j] = uc(j + 1);
    }
    
    double a1(0.0), a2(0.0), tmp;
    for (int n = 0; n < nfaces; ++n) {  
      int f = faces[n];
      const AmanziGeometry::Point& xf = mymesh->face_centroid(f);
      double area = mymesh->face_area(f);
       
      tmp = vf[n].Value(xf) - grad * (xf - xc);
      a1 += tmp * area;
      a2 += area;
    }

    uc(0) = a1 / a2;
  } else if (order_ >= 2) {
    v4 = integrals_.poly().coefs();
    basis.ChangeBasisMyToNatural(v4);
    v4.Reshape(nd);
    uc(0) = vdof(row) - (v4 * v5) / volume;
  }

  // calculate L2 projector
  if (type == ProjectorType::L2 && ndof_c > 0) {
    v5(0) = uc(0);

    DenseMatrix M, M2;
    DenseVector v6(nd - ndof_c);
    NumericalIntegration<AmanziMesh::Mesh> numi(mymesh);

    GrammMatrix(numi, order_, integrals_, basis, M);
    M2 = M.SubMatrix(ndof_c, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    const DenseVector& v3 = moments->coefs();
    for (int n = 0; n < ndof_c; ++n) {
      v4(n) = v3(n) * mymesh->cell_volume(c);
    }

    for (int n = 0; n < nd - ndof_c; ++n) {
      v4(ndof_c + n) = v6(n);
    }

    M.Inverse();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mymesh, c, order_, v5);
  }

  // set correct origin 
  uc.set_origin(xc);
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

