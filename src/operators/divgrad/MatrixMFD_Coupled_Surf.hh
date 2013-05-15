/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  Class to form the 2x n_lambda system for coupling energy and flow on both
  the surface and subsurface.

*/

#ifndef OPERATORS_MATRIX_COUPLED_MFD_SURF_HH_
#define OPERATORS_MATRIX_COUPLED_MFD_SURF_HH_


#include "MatrixMFD_TPFA.hh"
#include "MatrixMFD_Coupled.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Coupled_Surf : public MatrixMFD_Coupled {

 public:
  MatrixMFD_Coupled_Surf(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                       const Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh);

  MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other);

  virtual void ComputeSchurComplement(const Epetra_MultiVector& Ccc,
          const Epetra_MultiVector& Dcc,
          const Teuchos::Ptr<const Epetra_MultiVector>& Ccc_surf=Teuchos::null,
          const Teuchos::Ptr<const Epetra_MultiVector>& Dcc_surf=Teuchos::null);

  void SetSurfaceOperators(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                           const Teuchos::RCP<MatrixMFD_TPFA>& surface_B) {
    surface_A_ = surface_A;
    surface_B_ = surface_B;
  }

  void GetSurfaceOperators(Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                           Teuchos::RCP<MatrixMFD_TPFA>& surface_B) {
    surface_A = surface_A_;
    surface_B = surface_B_;
  }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_B_;

};


} // namespace
} // namespace


#endif
