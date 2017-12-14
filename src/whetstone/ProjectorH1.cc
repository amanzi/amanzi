/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Algorithms underpinning elliptic projectors.
*/

#include "MFD3D_CrouzeixRaviart.hh"
#include "NumericalIntegration.hh"
#include "ProjectorH1.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Enegy projector on space of linear polynomials in cell c.
* Uniquness require to specify its value at cell centroid.
****************************************************************** */
void ProjectorH1::HarmonicP1_Cell(
    int c, const AmanziGeometry::Point& p0,
    const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  double vol = mesh_->cell_volume(c);

  // create zero vector polynomial
  uc.resize(d_);
  for (int i = 0; i < d_; ++i) { 
    uc[i].Reshape(d_, 1, true);
  }

  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    for (int i = 0; i < d_; ++i) {
      double tmp = vf[n][i].Value(xf) * dirs[n] / vol;

      for (int j = 0; j < d_; ++j) {
        uc[i].monomials(1).coefs()[j] += tmp * normal[j];
      }
    }
  }

  // set projector zero term to the given p0
  const AmanziGeometry::Point& xc0 = mesh_->cell_centroid(c);
  AmanziGeometry::Point zero(d_);

  for (int i = 0; i < d_; ++i) {
    uc[i](0, 0) = p0[i];
    uc[i].set_origin(xc0);
    uc[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Energy projector on space of linear polynomials in face f.
* Uniquness require to specy its value at cell centroid.
****************************************************************** */
void ProjectorH1::HarmonicP1_Face(
    int f, const AmanziGeometry::Point& p0,
    const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const
{
  Entity_ID_List edges;
  std::vector<int> dirs;

  mesh_->face_get_edges_and_dirs(f, &edges, &dirs);
  int nedges = edges.size();

  double area = mesh_->face_area(f);
  AmanziGeometry::Point fnormal = mesh_->face_normal(f);
  fnormal /= norm(fnormal);

  // create zero vector polynomial
  uf.resize(d_);
  for (int i = 0; i < d_; ++i) { 
    uf[i].Reshape(d_, 1, true);
  }

  AmanziGeometry::Point enormal(d_);

  for (int n = 0; n < nedges; ++n) {  
    int e = edges[n];
    const AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
    const AmanziGeometry::Point& tau = mesh_->edge_vector(e);

    enormal = tau^fnormal;

    for (int i = 0; i < d_; ++i) {
      double tmp = ve[n][i].Value(xe) * dirs[n] / area;

      for (int j = 0; j < d_; ++j) {
        uf[i](1, j) += tmp * enormal[j];
      }
    }
  }

  // fix the constant value
  const AmanziGeometry::Point& xf0 = mesh_->face_centroid(f);
  AmanziGeometry::Point zero(d_);

  for (int i = 0; i < d_; ++i) {
    uf[i](0, 0) = p0[i];
    uf[i].set_origin(xf0);
    uf[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Enegy projector on space of polynomials of order k in cell c.
* Uniquness require to specify its value at cell centroid.
****************************************************************** */
void ProjectorH1::HarmonicPk_Cell(
    int c, const AmanziGeometry::Point& p0, int order,
    const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const
{
  ASSERT(d_ == 2);

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, R, Gpoly, A;
  MFD3D_CrouzeixRaviart mfd(mesh_);

  T(0, 0) = 1.0;
  // mfd.ModifyStabilityScalingFactor(2.0);
  mfd.StiffnessMatrixHO(c, order, T, R, Gpoly, A);  

  // number of degrees of freedom
  Polynomial pf(d_ - 1, order - 1);
  int nd = Gpoly.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  DenseMatrix Acf = A.SubMatrix(ndof_f, ndof, 0, ndof_f);
  DenseMatrix Acc = A.SubMatrix(ndof_f, ndof, ndof_f, ndof);
  Acc.Inverse();
  
  // create zero vector polynomial
  uc.resize(d_);
  for (int i = 0; i < d_; ++i) { 
    uc[i].Reshape(d_, order, true);
  }

  // calculate DOFs for boundary polynomial
  DenseVector vdof(ndof);
  std::vector<const Polynomial*> polys(2);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  NumericalIntegration numi(mesh_);

  for (int i = 0; i < d_; ++i) {
    int row(0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      if (d_ == 2) {
        tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
      }

      polys[0] = &(vf[n][i]);

      for (auto it = pf.begin(); it.end() <= pf.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 1, index);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys);
        row++;
      }
    }

    // harmonic extension inside cell
    DenseVector v1(ndof_f), v2(ndof_c), v3(ndof_c);

    for (int n = 0; n < ndof_f; ++n) {
      v1(n) = vdof(n);
    }

    Acf.Multiply(v1, v2, false);
    Acc.Multiply(v2, v3, false);

    for (int n = 0; n < ndof_c; ++n) {
      vdof(row + n) = -v3(n);
    }

    // calculate polynomial coefficients
    DenseVector v4(nd), v5(nd);
    R.Multiply(vdof, v4, true);
    Gpoly.Multiply(v4, v5, false);

    uc[i].set_coefs(v5);
  }

  // set projector zero term to the given p0
  const AmanziGeometry::Point& xc0 = mesh_->cell_centroid(c);
  AmanziGeometry::Point zero(d_);

  for (int i = 0; i < d_; ++i) {
    uc[i](0, 0) = p0[i];
    uc[i].set_origin(xc0);
    uc[i].ChangeOrigin(zero);
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

