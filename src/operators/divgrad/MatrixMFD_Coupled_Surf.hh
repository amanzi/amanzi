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
                         const Teuchos::RCP<const AmanziMesh::Mesh> mesh);
  MatrixMFD_Coupled_Surf(const MatrixMFD_Coupled_Surf& other);

  virtual void SetSubBlocks(const Teuchos::RCP<MatrixMFD>& blockA,
                            const Teuchos::RCP<MatrixMFD>& blockB) {
    MatrixMFD_Coupled::SetSubBlocks(blockA, blockB);
  }

  virtual void ComputeSchurComplement();

  virtual void SetOffDiagonals(const Teuchos::RCP<const Epetra_MultiVector>& Ccc,
          const Teuchos::RCP<const Epetra_MultiVector>& Dcc,
          const Teuchos::RCP<const Epetra_MultiVector>& Ccc_surf=Teuchos::null,
          const Teuchos::RCP<const Epetra_MultiVector>& Dcc_surf=Teuchos::null,
          double scaling=1.) {
    MatrixMFD_Coupled::SetOffDiagonals(Ccc,Dcc,scaling);
    if (Ccc_surf == Teuchos::null) {
      Teuchos::RCP<Epetra_MultiVector> Ccc_s = 
          Teuchos::rcp(new Epetra_MultiVector(surface_mesh_->cell_map(false),1));
      Ccc_s->PutScalar(0.);
      Ccc_surf_ = Ccc_s;
    } else {
      Ccc_surf_ = Ccc_surf;
    }

    if (Dcc_surf == Teuchos::null) {
      Teuchos::RCP<Epetra_MultiVector> Dcc_s = 
          Teuchos::rcp(new Epetra_MultiVector(surface_mesh_->cell_map(false),1));
      Dcc_s->PutScalar(0.);
      Dcc_surf_ = Dcc_s;
    } else {
      Dcc_surf_ = Dcc_surf;
    }
  }
  
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
  Teuchos::RCP<const AmanziMesh::Mesh> surface_mesh_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_A_;
  Teuchos::RCP<MatrixMFD_TPFA> surface_B_;

  Teuchos::RCP<const Epetra_MultiVector> Ccc_surf_;
  Teuchos::RCP<const Epetra_MultiVector> Dcc_surf_;

  Teuchos::RCP<const Epetra_Import> surf_importer_;
  Teuchos::RCP<const Epetra_Map> surf_map_in_subsurf_;
  
  bool dump_schur_;
};


} // namespace
} // namespace


#endif
