/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Bernardi Raugel element: velocity space vectors at nodes and 
  normal components on faces.
*/

#ifndef AMANZI_MFD3D_BERNARDI_RAUGEL_HH_
#define AMANZI_MFD3D_BERNARDI_RAUGEL_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_BernardiRaugel : public MFD3D {
 public:
  MFD3D_BernardiRaugel(const Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
    : MFD3D(mesh) {};

  // required methods
  // -- schemas
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- other matrices
  virtual int DivergenceMatrix(int c, DenseMatrix& A) override;
  virtual int AdvectionMatrix(int c, const std::vector<AmanziGeometry::Point>& u, DenseMatrix& A) override;

 private:
  static RegisteredFactory<MFD3D_BernardiRaugel> factory_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

