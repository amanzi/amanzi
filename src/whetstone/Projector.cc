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

#include "DG_Modal.hh"
#include "MFD3D_CrouzeixRaviart.hh"
#include "MFD3D_Lagrange.hh"
#include "NumericalIntegration.hh"
#include "Projector.hh"
#include "VectorPolynomial.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Energy projector on the space of linear polynomials in cell c.
****************************************************************** */
void Projector::HarmonicCell_CR1(
    int c, const std::vector<VectorPolynomial>& vf, VectorPolynomial& uc) const
{
  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double vol = mesh_->cell_volume(c);

  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, 1, true);
  }

  for (int n = 0; n < nfaces; ++n) {  
    int f = faces[n];
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);

    for (int i = 0; i < dim; ++i) {
      double tmp = vf[n][i].Value(xf) * dirs[n] / vol;

      for (int j = 0; j < d_; ++j) {
        uc[i].monomials(1).coefs()[j] += tmp * normal[j];
      }
    }
  }

  // calculate projector's low-order term
  AmanziGeometry::Point grad(d_), zero(d_);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < d_; ++j) {
      grad[j] = uc[i](1, j);
    }
    
    double a1(0.0), a2(0.0), tmp;
    for (int n = 0; n < nfaces; ++n) {  
      int f = faces[n];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      double area = mesh_->face_area(f);
       
      tmp = vf[n][i].Value(xf) - grad * (xf - xc);
      a1 += tmp * area;
      a2 += area;
    }

    uc[i](0, 0) = a1 / a2;
  }

  // clean-up: fix the origin
  for (int i = 0; i < dim; ++i) {
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Energy projector on space of linear polynomials in face f.
* Uniqueness requires to specify projector's value at face centroid.
****************************************************************** */
void Projector::HarmonicFace_CR1(
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
* Energy projector on space of polynomials of order k in cell c.
* Uniqueness require to specify its value at cell centroid.
****************************************************************** */
void Projector::GenericCell_CRk_(
    int c, int order, const std::vector<VectorPolynomial>& vf,
    const Projector::Type type, bool is_harmonic, 
    const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const
{
  ASSERT(d_ == 2);

  Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

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

  DenseMatrix Acf, Acc;
  if (ndof_c > 0 && is_harmonic) {
    Acf = A.SubMatrix(ndof_f, ndof, 0, ndof_f);
    Acc = A.SubMatrix(ndof_f, ndof, ndof_f, ndof);
    Acc.Inverse();
  }
  
  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, order, true);
  }

  // calculate DOFs for boundary polynomial
  DenseVector vdof(ndof);
  std::vector<const Polynomial*> polys(2);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);
  NumericalIntegration numi(mesh_);

  for (int i = 0; i < dim; ++i) {
    int row(0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      if (d_ == 2) {
        tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
      }

      polys[0] = &(vf[n][i]);

      for (auto it = pf.begin(); it.end() <= pf.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }

    // harmonic extension inside cell
    if (ndof_c > 0 && is_harmonic) {
      DenseVector v1(ndof_f), v2(ndof_c), v3(ndof_c);

      for (int n = 0; n < ndof_f; ++n) {
        v1(n) = vdof(n);
      }

      Acf.Multiply(v1, v2, false);
      Acc.Multiply(v2, v3, false);

      moments->Reshape(ndof_c);
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = -v3(n);
        (*moments)(n) = -v3(n);
      }
    }
    else if (ndof_c > 0 && !is_harmonic) {
      ASSERT(ndof_c == moments->NumRows());
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = (*moments)(n);
      }
    }

    // calculate polynomial coefficients
    DenseVector v4(nd), v5(nd);
    R.Multiply(vdof, v4, true);
    Gpoly.Multiply(v4, v5, false);

    uc[i].SetPolynomialCoefficients(v5);
    numi.ChangeBasisNaturalToRegular(c, uc[i]);

    // calculate the constant value
    if (order == 1) {
      AmanziGeometry::Point grad(d_), zero(d_);
      for (int j = 0; j < d_; ++j) {
        grad[j] = uc[i](1, j);
      }
    
      double a1(0.0), a2(0.0), tmp;
      for (int n = 0; n < nfaces; ++n) {  
        int f = faces[n];
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        double area = mesh_->face_area(f);
       
        tmp = vf[n][i].Value(xf) - grad * (xf - xc);
        a1 += tmp * area;
        a2 += area;
      }

      uc[i](0, 0) = a1 / a2;
    } else if (order >= 2) {
      mfd.integrals().GetPolynomialCoefficients(v4);
      v4.Reshape(nd);
      uc[i](0, 0) = vdof(row) - (v4 * v5) / volume;
    }

    // change origin from centroid to zero
    AmanziGeometry::Point zero(d_);
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Energy projector on space of polynomials of order k in cell c.
* Uniqueness require to specify its value at cell centroid.
****************************************************************** */
void Projector::GenericCell_Pk_(
    int c, int order, const std::vector<VectorPolynomial>& vf, 
    const Projector::Type type, bool is_harmonic,
    const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const
{
  ASSERT(d_ == 2);

  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  double volume = mesh_->cell_volume(c);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, R, Gpoly, A;
  MFD3D_Lagrange mfd(mesh_);

  T(0, 0) = 1.0;
  // mfd.ModifyStabilityScalingFactor(1.02);
  mfd.StiffnessMatrixHO(c, order, T, R, Gpoly, A);  

  // number of degrees of freedom
  Polynomial pf;
  if (order > 1) {
    pf.Reshape(d_ - 1, order - 2);
  }

  int nd = Gpoly.NumRows();
  int ndf = pf.size();
  int ndof = A.NumRows();

  int ndof_f(nnodes + nfaces * ndf);
  int ndof_c(ndof - ndof_f);

  // auxiliary data for optional calculation of internal moments
  DenseMatrix Acf, Acc;
  if (ndof_c > 0 && is_harmonic) {
    Acf = A.SubMatrix(ndof_f, ndof, 0, ndof_f);
    Acc = A.SubMatrix(ndof_f, ndof, ndof_f, ndof);
    Acc.Inverse();
  }
  
  // create zero vector polynomial
  int dim = vf[0].size();
  uc.resize(dim);
  for (int i = 0; i < dim; ++i) { 
    uc[i].Reshape(d_, order, true);
  }

  DenseVector vdof(ndof);
  std::vector<const Polynomial*> polys(2);
  NumericalIntegration numi(mesh_);

  AmanziGeometry::Point xv(d_);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);

  for (int i = 0; i < dim; ++i) {
    int row(nnodes);

    // calculate DOFs on boundary
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];

      Entity_ID_List face_nodes;
      mesh_->face_get_nodes(f, &face_nodes);
      int nfnodes = face_nodes.size();

      for (int j = 0; j < nfnodes; j++) {
        int v = face_nodes[j];
        mesh_->node_get_coordinates(v, &xv);

        int pos = WhetStone::FindPosition(v, nodes);
        vdof(pos) = vf[n][i].Value(xv);
      }

      if (order > 1) { 
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
        double area = mesh_->face_area(f);

        // local coordinate system with origin at face centroid
        const AmanziGeometry::Point& normal = mesh_->face_normal(f);
        if (d_ == 2) {
          tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
        }

        polys[0] = &(vf[n][i]);

        for (auto it = pf.begin(); it.end() <= pf.end(); ++it) {
          const int* index = it.multi_index();
          Polynomial fmono(d_ - 1, index, 1.0);
          fmono.InverseChangeCoordinates(xf, tau);  

          polys[1] = &fmono;

          vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
          row++;
        }
      }
    }

    // calculate DOFs inside cell using
    // -- either harmonic extension inside cell
    if (ndof_c > 0 && is_harmonic) {
      DenseVector v1(ndof_f), v2(ndof_c), v3(ndof_c);

      for (int n = 0; n < ndof_f; ++n) {
        v1(n) = vdof(n);
      }

      Acf.Multiply(v1, v2, false);
      Acc.Multiply(v2, v3, false);

      moments->Reshape(ndof_c);
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = -v3(n);
        (*moments)(n) = -v3(n);
      }
    }
    // -- or copy moments from input data
    else if (ndof_c > 0 && !is_harmonic) {
      ASSERT(ndof_c == moments->NumRows());
      for (int n = 0; n < ndof_c; ++n) {
        vdof(row + n) = (*moments)(n);
      }
    }

    // calculate polynomial coefficients (in vector v5)
    DenseVector v4(nd), v5(nd);
    R.Multiply(vdof, v4, true);
    Gpoly.Multiply(v4, v5, false);

    uc[i].SetPolynomialCoefficients(v5);
    numi.ChangeBasisNaturalToRegular(c, uc[i]);

    // calculate the constant value for elliptic projector
    if (order == 1) {
      AmanziGeometry::Point grad(d_), zero(d_);
      for (int j = 0; j < d_; ++j) {
        grad[j] = uc[i](1, j);
      }
    
      double a1(0.0), a2(0.0), tmp;
      for (int n = 0; n < nfaces; ++n) {  
        int f = faces[n];
        const AmanziGeometry::Point& xf = mesh_->face_centroid(f);
        double area = mesh_->face_area(f);
       
        tmp = vf[n][i].Value(xf) - grad * (xf - xc);
        a1 += tmp * area;
        a2 += area;
      }

      uc[i](0, 0) = a1 / a2;
    } else if (order >= 2) {
      mfd.integrals().poly().GetPolynomialCoefficients(v4);
      v4.Reshape(nd);
      uc[i](0, 0) = vdof(row) - (v4 * v5) / volume;
    }

    // calculate L2 projector
    if (type == Type::L2 && ndof_c > 0) {
      v5(0) = uc[i](0, 0);

      DG_Modal dg(order, mesh_);
      dg.set_basis(TAYLOR_BASIS_NATURAL);

      DenseMatrix M, M2;
      DenseVector v6(nd - ndof_c);
      dg.MassMatrix(c, T, mfd.integrals(), M);

      M2 = M.SubMatrix(ndof_c, nd, 0, nd);
      M2.Multiply(v5, v6, false);

      for (int n = 0; n < ndof_c; ++n) {
        v4(n) = (*moments)(n) * mesh_->cell_volume(c);
      }

      for (int n = 0; n < nd - ndof_c; ++n) {
        v4(ndof_c + n) = v6(n);
      }

      M.Inverse();
      M.Multiply(v4, v5, false);

      uc[i].SetPolynomialCoefficients(v5);
      numi.ChangeBasisNaturalToRegular(c, uc[i]);
    }

    // set origin to zero
    AmanziGeometry::Point zero(d_);
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* L2 projector form serendipidy space.
****************************************************************** */
void Projector::L2Cell_SerendipityPk(
    int c, int order, const std::vector<VectorPolynomial>& vf,
    const std::shared_ptr<DenseVector>& moments, VectorPolynomial& uc) const
{
  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, R, Gpoly, A;
  MFD3D_Lagrange mfd(mesh_);

  T(0, 0) = 1.0;
  mfd.H1consistencyHO(c, order, T, N, R, Gpoly, A);  

  NumericalIntegration numi(mesh_);

  // number of degrees of freedom
  Polynomial pc;
  if (order > 1) {
    pc.Reshape(d_, order - 2);
  }

  int nd = Gpoly.NumRows();
  int ndof = A.NumRows();
  int ndof_c(pc.size());
  int ndof_f(ndof - ndof_c);

  // extract submatrix
  DenseMatrix Ns, NN(nd, nd);
  Ns = N.SubMatrix(0, ndof_f, 0, nd);

  NN.Multiply(Ns, Ns, true);
  NN.Inverse();

  // calculate degrees of freedom
  const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
  DenseVector vdof(ndof_f), v1(nd), v2(nd);

  int dim = vf[0].size();
  uc.resize(dim);

  for (int i = 0; i < dim; ++i) {
    CalculateDOFsOnBoundary_Pk_(c, order, vf, vdof, i);

    Ns.Multiply(vdof, v1, true);
    NN.Multiply(v1, v2, false);

    uc[i].Reshape(d_, order, true);
    uc[i].SetPolynomialCoefficients(v2);
    numi.ChangeBasisNaturalToRegular(c, uc[i]);

    // set origin to zero
    AmanziGeometry::Point zero(d_);
    uc[i].set_origin(xc);
    uc[i].ChangeOrigin(zero);
  }
}


/* ******************************************************************
* Calculate degrees of freedom in 2D.
****************************************************************** */
void Projector::CalculateDOFsOnBoundary_Pk_(
    int c, int order, const std::vector<VectorPolynomial>& vf,
    DenseVector& vdof, int i) const
{
  Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  mesh_->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  std::vector<const Polynomial*> polys(2);
  NumericalIntegration numi(mesh_);

  AmanziGeometry::Point xv(d_);
  std::vector<AmanziGeometry::Point> tau(d_ - 1);

  // number of moments of faces
  Polynomial pf;
  if (order > 1) {
    pf.Reshape(d_ - 1, order - 2);
  }

  int row(nnodes);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    Entity_ID_List face_nodes;
    mesh_->face_get_nodes(f, &face_nodes);
    int nfnodes = face_nodes.size();

    for (int j = 0; j < nfnodes; j++) {
      int v = face_nodes[j];
      mesh_->node_get_coordinates(v, &xv);

      int pos = WhetStone::FindPosition(v, nodes);
      vdof(pos) = vf[n][i].Value(xv);
    }

    if (order > 1) { 
      const AmanziGeometry::Point& xf = mesh_->face_centroid(f); 
      double area = mesh_->face_area(f);

      // local coordinate system with origin at face centroid
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      if (d_ == 2) {
        tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
      }

      polys[0] = &(vf[n][i]);

      for (auto it = pf.begin(); it.end() <= pf.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d_ - 1, index, 1.0);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

