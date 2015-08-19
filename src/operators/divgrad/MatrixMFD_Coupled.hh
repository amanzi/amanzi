/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Coupled takes two MatrixMFD objects, along with the cell coupling
  terms, and forms a coupled system that is 2x cell + 2x face sqare.

  MatrixMFD provides for a system:

      [  Acc   Acf  ]
      [  Afc   Aff  ]

  where Acc is size ncells x ncells and Aff is size nfaces x nfaces.  The key
  to solving this system is that Acc is diagonal, allowing us to precondition
  with the Schur complement:

      [  Acc   Acf  ]
      [   0    Sff  ],   Sff = Aff - Afc * Acc^-1 * Acf

  For many of the systems needed in ATS, we have two, coupled equations whose
  preconditioner looks something like:

      [  Acc   Acf    Ccc    0  ]
      [  Afc   Aff     0     0  ]
      [  Dcc    0     Bcc   Bcf ]
      [   0     0     Bfc   Bff ]

  where Acc, Bcc, Ccc, and Dcc are all diagonal.  Effectively these four
  diagonal matrices form a block-diagonal matrix where each cell corresponds
  to a 2x2 dense matrix on the diagonal.  The same Schur complement trick can
  then be applied:

     S_(2f x 2f) = [ Aff   0  ]  _  [ Afc  0  ][ Acc  Ccc ]^-1 [ Acf  0  ]
                   [  0   Bff ]     [  0  Bfc ][ Dcc  Bcc ]    [  0  Bcf ]

  This class forms this matrix.
 */


#ifndef OPERATORS_MATRIX_COUPLED_MFD_HH_
#define OPERATORS_MATRIX_COUPLED_MFD_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_FEVbrMatrix.h"

#include "TreeVector.hh"
#include "TreeMatrix.hh"
#include "Preconditioner.hh"

#include "MatrixMFD.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Coupled : public TreeMatrix {

 public:
  MatrixMFD_Coupled(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  MatrixMFD_Coupled(const MatrixMFD_Coupled& other);

  void InitializeFromPList_();

  virtual void SetSubBlocks(const Teuchos::RCP<MatrixMFD>& blockA,
                            const Teuchos::RCP<MatrixMFD>& blockB);

  virtual void SetAdvectiveBlock(const Teuchos::RCP<MatrixMFD>& adv_block) {
    adv_block_ = adv_block;
  }

  virtual void SetOffDiagonals(const Teuchos::RCP<const Epetra_MultiVector>& Ccc,
                       const Teuchos::RCP<const Epetra_MultiVector>& Dcc,
                       double scaling) {
    scaling_ = scaling;
    Ccc_ = Ccc;
    Dcc_ = Dcc;
    MarkLocalMatricesAsChanged_();
  }

  Teuchos::RCP<const Epetra_FEVbrMatrix> Schur() {
    return P2f2f_;
  }
  Teuchos::RCP<const Epetra_FEVbrMatrix> Aff() {
    return A2f2f_;
  }

  // TreeMatrix stuff FIX ME!
  virtual const TreeVectorSpace& DomainMap() const {
    return *space_; }

  // Vector space of the Matrix's range.
  virtual const TreeVectorSpace& RangeMap() const {
    return *space_; }

  // Virtual copy constructor.
  virtual Teuchos::RCP<TreeMatrix> Clone() const {
    return Teuchos::rcp(new MatrixMFD_Coupled(*this));
  }

  virtual int Apply(const TreeVector& X,
                     TreeVector& Y) const;
  virtual int ApplyInverse(const TreeVector& X,
                            TreeVector& Y) const;

  virtual void Apply(const TreeVector& X,
                     const Teuchos::Ptr<TreeVector>& Y) const {
    Apply(X,*Y); }
  virtual void ApplyInverse(const TreeVector& X,
                            const Teuchos::Ptr<TreeVector>& Y) const {
    ApplyInverse(X, *Y); }

  virtual void SymbolicAssembleGlobalMatrices();
  virtual void InitPreconditioner() {}
  virtual void UpdateConsistentFaceCorrection(const TreeVector& u,
          const Teuchos::Ptr<TreeVector>& Pu);

 protected:
  virtual void MarkLocalMatricesAsChanged_() {
    assembled_schur_ = false;
    assembled_operator_ = false;
  }

  virtual void AssembleSchur_() const;
  virtual void AssembleAff_() const;
  virtual void UpdatePreconditioner_() const;

 protected:
  // mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;

  // domain/range space
  Teuchos::RCP<TreeVectorSpace> space_;

  // sub-blocks
  Teuchos::RCP<MatrixMFD> blockA_;
  Teuchos::RCP<MatrixMFD> blockB_;
  Teuchos::RCP<MatrixMFD> adv_block_;

  // off diagonal blocks
  Teuchos::RCP<const Epetra_MultiVector> Ccc_;
  Teuchos::RCP<const Epetra_MultiVector> Dcc_;
  double scaling_;

  // local matrices
  mutable std::vector<Teuchos::SerialDenseMatrix<int, double> > A2c2c_cells_Inv_;

  // global matrices
  mutable Teuchos::RCP<Epetra_VbrMatrix> A2f2c_;
  mutable Teuchos::RCP<Epetra_VbrMatrix> A2c2f_;
  mutable Teuchos::RCP<Epetra_FEVbrMatrix> P2f2f_;
  mutable Teuchos::RCP<Epetra_FEVbrMatrix> A2f2f_;

  // maps
  Teuchos::RCP<const Epetra_BlockMap> double_fmap_;
  Teuchos::RCP<const Epetra_BlockMap> double_cmap_;
  Teuchos::RCP<const Epetra_BlockMap> double_fmap_wghost_;

  // flags
  mutable bool assembled_operator_;
  mutable bool assembled_schur_;
  mutable bool is_schur_created_;
  mutable bool is_operator_created_;
  bool dump_schur_;

  // preconditioner for Schur complement
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> S_pc_;

  // verbose object
  Teuchos::RCP<VerboseObject> vo_;

};

} // namespace
} // namespace

#endif
