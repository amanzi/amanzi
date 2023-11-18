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

  The mimetic finite difference method for elasticity with weak symmetry.
  Stress: linear normal component per face (d*d dofs)
    (1) order is: zero moments, d first moments
  Displacement: 1 value per cell (d dofs)
  Rotations:    1 value per cell (d dofs)
*/

#ifndef AMANZI_MFD3D_ELASTICITY_WEAK_SYMMETRY_HH_
#define AMANZI_MFD3D_ELASTICITY_WEAK_SYMMETRY_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_ElasticityWeakSymmetry : public MFD3D {
 public:
  MFD3D_ElasticityWeakSymmetry(const Teuchos::ParameterList& plist,
                               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh){};

  // required methods
  // -- schema for stress, displacement and rotations
  virtual std::vector<SchemaItem> schema() const override
  {
    std::vector<SchemaItem> items;
    items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, d_));
    items.push_back(std::make_tuple(AmanziMesh::CELL, DOF_Type::SCALAR, d_));
    items.push_back(std::make_tuple(AmanziMesh::NODE, DOF_Type::SCALAR, d_ * (d_ - 1) / 2));
    return items;
  }

  // -- mass matrices
  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override;

  // -- stiffness matrix
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- divergence matrix
  virtual int DivergenceMatrix(int c, DenseMatrix& B) override;

  // rotation metrix
  void RotationMatrix(int c, DenseMatrix& G);

 private:
  void set_index_(int d, int k, int* index);

  void TripleMatrixProduct_(const DenseVector& ML,
                            const DenseMatrix& Minv,
                            const DenseVector& MR,
                            DenseMatrix& A,
                            int i0,
                            int j0);

  void TripleMatrixProduct_(const DenseVector& ML,
                            const DenseMatrix& Minv,
                            const DenseMatrix& MR,
                            DenseMatrix& A,
                            int i0,
                            int j0);

  void TripleMatrixProduct_(const DenseMatrix& ML,
                            const DenseMatrix& Minv,
                            const DenseMatrix& MR,
                            DenseMatrix& A,
                            int i0,
                            int j0);

  static RegisteredFactory<MFD3D_ElasticityWeakSymmetry> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
