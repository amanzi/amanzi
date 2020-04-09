/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_WHETSTONE_MFD3D_DIFFUSION_EDGE_HH_
#define AMANZI_WHETSTONE_MFD3D_DIFFUSION_EDGE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "Tensor.hh"
#include "MFD3D.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_Diffusion_Edge : public MFD3D {
 public:
  MFD3D_Diffusion_Edge(const Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh), InnerProduct(mesh){};
  ~MFD3D_Diffusion_Edge(){};

  // main methods
  // -- mass matrices (not implemented)
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N,
                            DenseMatrix& Mc, bool symmetry)
  {
    return -1;
  }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) { return -1; }

  // -- inverse mass matrix (not implemented)
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R,
                                   DenseMatrix& Wc, bool symmetry)
  {
    return -1;
  }
  virtual int MassMatrixInverse(int c, const Tensor& K, DenseMatrix& W)
  {
    return -1;
  }

  // -- stiffness matrix
  virtual int
  H1consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A);

 private:
  static RegisteredFactory<MFD3D_Diffusion_Edge> factory_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
