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
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviart : public MFD3D {
 public:
  MFD3D_CrouzeixRaviart(const Teuchos::ParameterList& plist,
                        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override
  {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::Entity_kind::FACE, DOF_Type::SCALAR, 1));
  }

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(mesh_, c, ve, vf, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& uc) override
  {
    ProjectorCell_(mesh_, c, ve, vf, uc);
  }

  virtual void H1Face(int f,
                      const std::vector<Polynomial>& ve,
                      const Polynomial* moments,
                      Polynomial& vf) override;

  // access / setup
  // -- integrals of monomials in high-order schemes could be reused
  const DenseMatrix& G() const { return G_; }
  const DenseMatrix& R() const { return R_; }

 private:
  // efficient implementation of low-order elliptic projectors
  void ProjectorCell_(const Teuchos::RCP<const AmanziMesh::Mesh>& mymesh,
                      int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      Polynomial& uc);

 protected:
  DenseMatrix R_, G_;

 private:
  static RegisteredFactory<MFD3D_CrouzeixRaviart> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
