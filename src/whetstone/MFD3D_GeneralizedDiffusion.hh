/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mimetic schemes for generalized polyhedra.
*/

#ifndef AMANZI_MFD3D_GENERALIZED_DIFFUSION_HH_
#define AMANZI_MFD3D_GENERALIZED_DIFFUSION_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_GeneralizedDiffusion : public MFD3D { 
 public:
  MFD3D_GeneralizedDiffusion(const Teuchos::ParameterList& plist,
                             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : BilinearForm(mesh) {};

  // required member functions
  // -- schema for this element
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrices
  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override; 
  int MassMatrixOptimized(int c, const Tensor& T, DenseMatrix& M); 

  int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) override; 

  // -- stiffness matrices
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // -- advection matrices
  virtual int DivergenceMatrix(int c, DenseMatrix& A) override;

 private:
  void CurvedFaceGeometry_(int f, int dirs,
                           std::vector<AmanziGeometry::Point>& vv, 
                           std::vector<AmanziGeometry::Point>& xm);

 private:
  static RegisteredFactory<MFD3D_GeneralizedDiffusion> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

