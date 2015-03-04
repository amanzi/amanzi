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
#include "EpetraMatrix.hh"
#include "EpetraMatrixDefault.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Coupled_Surf : public MatrixMFD_Coupled {

 public:
  MatrixMFD_Coupled_Surf(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh> mesh);
  MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other);

  virtual void SymbolicAssembleGlobalMatrices();

  virtual void SetOffDiagonals(const Teuchos::RCP<const Epetra_MultiVector>& Ccc,
          const Teuchos::RCP<const Epetra_MultiVector>& Dcc,
          const Teuchos::RCP<const Epetra_MultiVector>& Ccc_surf,
          const Teuchos::RCP<const Epetra_MultiVector>& Dcc_surf,
          double scaling);
  
  void SetSurfaceOperators(const Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                           const Teuchos::RCP<MatrixMFD_TPFA>& surface_B);

  void GetSurfaceOperators(Teuchos::RCP<MatrixMFD_TPFA>& surface_A,
                           Teuchos::RCP<MatrixMFD_TPFA>& surface_B) {
    surface_A = surface_A_;
    surface_B = surface_B_;
  }

  virtual int Apply(const TreeVector& X,
                    TreeVector& Y) const;
  virtual void Apply(const TreeVector& X,
                     const Teuchos::Ptr<TreeVector>& Y) const {
    Apply(X,*Y); }

  virtual void UpdateConsistentFaceCorrection(const TreeVector& u,
          const Teuchos::Ptr<TreeVector>& Pu);

 protected:
  virtual void AssembleAff_() const;
  virtual void AssembleSchur_() const;



 protected:
  bool dump_schur_; // this intentionally shadows the same name var in
		    // MatrixMFD_Coupled... you cannot both dump that
		    // one AND this one due to TRILINOS FAIL

  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_B_;

  Teuchos::RCP<EpetraMatrixDefault<Epetra_FEVbrMatrix> > A2f2f_op_;
  Teuchos::RCP<EpetraMatrix> A2f2f_solver_;
  
  Teuchos::RCP<const Epetra_MultiVector> Ccc_surf_;
  Teuchos::RCP<const Epetra_MultiVector> Dcc_surf_;
};


} // namespace
} // namespace


#endif
