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

#ifndef AMANZI_MFD3D_CROUZEIX_RAVIART_HH_
#define AMANZI_MFD3D_CROUZEIX_RAVIART_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MatrixPolynomial.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviart : public MFD3D {
 public:
  MFD3D_CrouzeixRaviart(const Teuchos::ParameterList& plist,
                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_CrouzeixRaviart(){};

  // required methods
  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N,
                            DenseMatrix& Mc, bool symmetry) override
  {
    Errors::Message msg(
      "L2 consistency is not implemented for Crouzeix-Raviart space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override
  {
    Errors::Message msg(
      "MassMatrix is not supported for Crouzeix-Raviart space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N,
                            DenseMatrix& Ac) override
  {
    if (order_ == 1 && !use_always_ho_) {
      return H1consistencyLO_(c, T, N, Ac);
    } else {
      return H1consistencyHO_(c, T, N, Ac);
    }
  }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override
  {
    if (order_ == 1 && !use_always_ho_)
      return StiffnessMatrixLO_(c, T, A);
    else
      return StiffnessMatrixHO_(c, T, A);
  }

  // -- projectors: base L2 and H1 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override
  {
    ProjectorCell_HO_(c, vf, ProjectorType::L2, moments, uc);
  }

  virtual void H1Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override
  {
    if (order_ == 1 && !use_always_ho_)
      ProjectorCell_LO_(c, vf, uc);
    else
      ProjectorCell_HO_(c, vf, ProjectorType::H1, moments, uc);
  }

  // additional miscaleneous projectors
  void L2GradientCell(int c, const std::vector<VectorPolynomial>& vf,
                      const std::shared_ptr<DenseVector>& moments,
                      MatrixPolynomial& uc)
  {
    ProjectorGradientCell_(c, vf, ProjectorType::L2, moments, uc);
  }

  void
  H1Face(int f, const AmanziGeometry::Point& p0,
         const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const;

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const PolynomialOnMesh& integrals() const { return integrals_; }
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }
  // -- modify internal parameters
  void set_use_always_ho(bool flag) { use_always_ho_ = flag; }

 private:
  // efficient implementation of low-order methods
  int H1consistencyLO_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int StiffnessMatrixLO_(int c, const Tensor& T, DenseMatrix& A);

  // high-order methods
  int H1consistencyHO_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  int StiffnessMatrixHO_(int c, const Tensor& T, DenseMatrix& A);

  // efficient implementation of low-order elliptic projectors
  void
  ProjectorCell_LO_(int c, const std::vector<Polynomial>& vf, Polynomial& uc);

  // generic code for multiple projectors
  void ProjectorCell_HO_(int c, const std::vector<Polynomial>& vf,
                         const ProjectorType type, const Polynomial* moments,
                         Polynomial& uc);

  void ProjectorGradientCell_(int c, const std::vector<VectorPolynomial>& vf,
                              const ProjectorType type,
                              const std::shared_ptr<DenseVector>& moments,
                              MatrixPolynomial& uc);

  // supporting routines
  void CalculateFaceDOFs_(int f, const Polynomial& vf, const Polynomial& pf,
                          DenseVector& vdof, int& row);

 protected:
  PolynomialOnMesh integrals_;
  DenseMatrix R_, G_;

 private:
  bool use_always_ho_;

  static RegisteredFactory<MFD3D_CrouzeixRaviart> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
