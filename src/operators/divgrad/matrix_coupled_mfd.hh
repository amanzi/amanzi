/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  MatrixCoupledMFD takes two MatrixMFD objects, along with the cell coupling
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

#include "ml_MultiLevelPreconditioner.h"
#include "Ifpack.h"
#include "Ifpack_ILU.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Teuchos_LAPACK.hpp"

#include "tree_vector.hh"
#include "matrix_mfd.hh"

namespace Amanzi {
namespace Operators {

class MatrixCoupledMFD : public Matrix {

 public:
  MatrixCoupledMFD(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<const AmanziMesh::Mesh> mesh);

  MatrixCoupledMFD(const MatrixCoupledMFD& other);

  void InitializeFromPList_();

  void SetSubBlocks(const Teuchos::RCP<MatrixMFD>& blockA,
                    const Teuchos::RCP<MatrixMFD>& blockB) {
    blockA_ = blockA;
    blockB_ = blockB;
  }

  virtual void Apply(const TreeVector& X,
                     const Teuchos::Ptr<TreeVector>& Y) const;
  virtual void ApplyInverse(const TreeVector& X,
                            const Teuchos::Ptr<TreeVector>& Y) const;

  virtual void ComputeSchurComplement(const Epetra_MultiVector& Ccc,
          const Epetra_MultiVector& Dcc);
  virtual void SymbolicAssembleGlobalMatrices();
  virtual void InitPreconditioner();
  virtual void UpdatePreconditioner();

 protected:
  // mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;

  // sub-blocks
  Teuchos::RCP<MatrixMFD> blockA_;
  Teuchos::RCP<MatrixMFD> blockB_;

  // local matrices
  std::vector<Teuchos::SerialDenseMatrix<int, double> > Aff_cells_;
  std::vector<Teuchos::SerialDenseMatrix<int, double> > A2c2c_cells_Inv_;

  // global matrices
  Teuchos::RCP<Epetra_VbrMatrix> A2f2c_;
  Teuchos::RCP<Epetra_FEVbrMatrix> P2f2f_;

  // maps
  Teuchos::RCP<const Epetra_BlockMap> double_fmap_;
  Teuchos::RCP<const Epetra_BlockMap> double_cmap_;
  Teuchos::RCP<const Epetra_BlockMap> double_fmap_wghost_;

  // flags
  bool is_matrix_constructed_;
  bool decoupled_;

  // preconditioning (This should be moved to Matrix?)
  enum PrecMethod { PREC_METHOD_NULL,
                    TRILINOS_ML,
                    TRILINOS_ILU,
                    TRILINOS_BLOCK_ILU,
                    HYPRE_AMG,
                    HYPRE_EUCLID,
                    HYPRE_PARASAILS };
  PrecMethod prec_method_;

  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec_;
  Teuchos::ParameterList ml_plist_;
  Teuchos::ParameterList coupled_plist_;

#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_Sff_;
  Teuchos::ParameterList hypre_plist_;
  int hypre_ncycles_, hypre_nsmooth_;
  double hypre_tol_, hypre_strong_threshold_;
#endif

  Teuchos::RCP<Ifpack_ILU> ilu_prec_;
  Teuchos::ParameterList ilu_plist_;

  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
  Teuchos::ParameterList ifp_plist_;


};

} // namespace
} // namespace

#endif
