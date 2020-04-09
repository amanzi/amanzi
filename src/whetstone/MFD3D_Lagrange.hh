/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Lagrange element: degrees of freedom are nodal values.
*/

#ifndef AMANZI_MFD3D_LAGRANGE_HH_
#define AMANZI_MFD3D_LAGRANGE_HH_

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Polynomial.hh"
//#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Lagrange : public MFD3D { 
 public:
  MFD3D_Lagrange(const Teuchos::ParameterList& plist,
                 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_Lagrange() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::NODE, DOF_Type::SCALAR, 1));
  }

  // -- mass matrices
  virtual int L2consistency(int c, const Tensor<>& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override {
    Errors::Message msg("L2 consistency is not implemented for Lagrange element.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }
  virtual int MassMatrix(int c, const Tensor<>& T, DenseMatrix& M) override {
    Errors::Message msg("Mass matrix is not implemented for Lagrange element.");
    Exceptions::amanzi_throw(msg);
    return 0;
  }

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor<>& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor<>& T, DenseMatrix& A) override;

  // -- projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, ve, vf, uc);
  }

  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_(c, ve, vf, uc);
  }

 private:
  void ProjectorCell_(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf, Polynomial& uc);

 private:
  static RegisteredFactory<MFD3D_Lagrange> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
