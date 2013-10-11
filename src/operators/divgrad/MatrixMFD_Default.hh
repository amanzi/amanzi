/*
  License: BSD
  Authors: Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov) (ATS version)
  MatrixMFD provides a mimetic discretization for the elliptic operator div K grad u.

*/

#ifndef OPERATORS_MATRIX_MFD_DEFAULT_HH_
#define OPERATORS_MATRIX_MFD_DEFAULT_HH_

#include "MatrixMFD_Factory.hh"
#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Default : public MatrixMFD {

 public:
  MatrixMFD_Default(Teuchos::ParameterList& plist,
		    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    MatrixMFD(plist, mesh) {}
    

private:
  // factory registration
  static RegisteredMatrixMFD_Factory<MatrixMFD_Default> reg_;
};

} // namespace
} // namespace

#endif
