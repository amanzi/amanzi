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

  The package uses the formula M = Mc + Ms, where matrix Mc is build from a
  consistency condition (Mc N = R) and matrix Ms is build from a stability
  condition  and Ms N = 0, to generate mass and stiffness matrices.
  The material properties are imbedded into the the matrix Mc.
*/

#ifndef AMANZI_MFD3D_ELASTICITY_HH_
#define AMANZI_MFD3D_ELASTICITY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class MFD3D_Elasticity : public MFD3D {
 public:
  MFD3D_Elasticity(const Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh) {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override
  {
    return std::vector<SchemaItem>(
      1, std::make_tuple(AmanziMesh::Entity_kind::NODE, DOF_Type::POINT, d_));
  }

  // -- stiffness matrices
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // optimization methods (mainly for research, since the maximum principle does not exists)
  int StiffnessMatrixOptimized(int c, const Tensor& T, DenseMatrix& A);
  int StiffnessMatrixMMatrix(int c, const Tensor& T, DenseMatrix& A);

  // projectors
  virtual void H1Cell(int c, const DenseVector& dofs, Tensor& vc) override;

 private:
  void MatrixMatrixProduct_(const DenseMatrix& A,
                            const DenseMatrix& B,
                            bool transposeB,
                            DenseMatrix& AB);

 private:
  DenseMatrix coefM_, R_;

  static RegisteredFactory<MFD3D_Elasticity> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
