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
  ~MFD3D_CrouzeixRaviart() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, 1));
  }

  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override {
    Errors::Message msg("L2 consistency is not implemented for Crouzeix-Raviart space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override {
    Errors::Message msg("MassMatrix is not supported for Crouzeix-Raviart space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- projectors: base L2 and H1 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, uc);
  }

  virtual void H1Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, vf, uc);
  }

  // additional miscaleneous projectors
  void H1Face(int f, const AmanziGeometry::Point& p0,
              const std::vector<VectorPolynomial>& ve, VectorPolynomial& uf) const;

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  // efficient implementation of low-order elliptic projectors
  void ProjectorCell_(int c, const std::vector<Polynomial>& vf, Polynomial& uc);

 protected:
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_CrouzeixRaviart> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

