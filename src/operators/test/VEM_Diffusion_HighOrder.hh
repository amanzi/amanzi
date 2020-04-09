/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Virtual element high-order conformal discretization of elliptic operator
  via a static condensation.
*/

#ifndef AMANZI_WHETSTONE_VEM3D_DIFFUSION_HIGH_ORDER_HH_
#define AMANZI_WHETSTONE_VEM3D_DIFFUSION_HIGH_ORDER_HH_

#include <tuple>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class VEM_Diffusion_HighOrder : public MFD3D {
 public:
  VEM_Diffusion_HighOrder(const Teuchos::ParameterList& plist,
                          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh), plist_(plist) {};
  ~VEM_Diffusion_HighOrder() {};

  // main methods 
  // -- symmetric schema
  virtual std::vector<SchemaItem> schema() const override {
    std::vector<SchemaItem> items;
    items.push_back(std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, 3));
    items.push_back(std::make_tuple(AmanziMesh::CELL, DOF_Type::SCALAR, 1));
    return items;
  }

  // -- mass matrices (not implemented)
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override { return -1; }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override { return -1; } 

  // -- inverse mass matrix (not implemented)
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) override { return -1; }
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W) override { return -1; } 

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac) override { return -1; }
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // post-processing method
  void UpdateFlux(int c, const Tensor& T, double** lambda, const double* p, double** flux);

 private:
  Teuchos::ParameterList plist_;

 private:
  static RegisteredFactory<VEM_Diffusion_HighOrder> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

