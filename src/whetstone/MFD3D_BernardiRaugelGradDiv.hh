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
*/

#ifndef AMANZI_MFD3D_BERNARDI_RAUGEL_GRAD_DIV_HH_
#define AMANZI_MFD3D_BERNARDI_RAUGEL_GRAD_DIV_HH_

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

class MFD3D_BernardiRaugelGradDiv : public MFD3D {
 public:
  MFD3D_BernardiRaugelGradDiv(const Teuchos::ParameterList& plist,
                              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : MFD3D(mesh) {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrices
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

 private:
  static RegisteredFactory<MFD3D_BernardiRaugelGradDiv> reg_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
