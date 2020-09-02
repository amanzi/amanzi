/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are mean values on faces.
*/

#ifndef AMANZI_MFD3D_CROUZEIX_RAVIART_HH_
#define AMANZI_MFD3D_CROUZEIX_RAVIART_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "SurfaceMiniMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviart : public MFD3D { 
 public:
  MFD3D_CrouzeixRaviart(const Teuchos::ParameterList& plist,
                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, 1));
  }

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, uc);
  }

  virtual void H1Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& vf) override ;

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  // efficient implementation of low-order elliptic projectors
  template<class MyMesh>
  void ProjectorCell_(const Teuchos::RCP<const MyMesh>& mymesh, 
                      int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf, Polynomial& uc);

 protected:
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_CrouzeixRaviart> factory_;
};


/* ******************************************************************
* Energy projector on the space of linear polynomials in cell c.
****************************************************************** */
template<class MyMesh>
void MFD3D_CrouzeixRaviart::ProjectorCell_(
    const Teuchos::RCP<const MyMesh>& mymesh, 
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf, Polynomial& uc)
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mymesh->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mymesh->cell_centroid(c);
  double vol = mymesh->cell_volume(c);

  // create zero vector polynomial
  uc.Reshape(d_, 1, true);

  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mymesh->face_centroid(f);
    const AmanziGeometry::Point& normal = mymesh->face_normal(f);

    double tmp = vf[n].Value(xf) * dirs[n] / vol;

    for (int j = 0; j < d_; ++j) {
      uc(1, j) += tmp * normal[j];
    }
  }

  // calculate projector's low-order term
  AmanziGeometry::Point grad(d_);
  for (int j = 0; j < d_; ++j) {
    grad[j] = uc(1, j);
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

  // set the correct origin
  uc.set_origin(xc);
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

