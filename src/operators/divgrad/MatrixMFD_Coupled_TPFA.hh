/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

  MatrixMFD_Coupled_TPFA takes two MatrixMFD_TPFA objects, along with
  the cell coupling terms, and forms a coupled system that is 2x cell
  + 2x face sqare.
 */


#ifndef OPERATORS_MATRIX_COUPLED_TPFA_HH_
#define OPERATORS_MATRIX_COUPLED_TPFA_HH_

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

#include "MatrixMFD_TPFA.hh"
#include "MatrixMFD_Coupled.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD_Coupled_TPFA : public MatrixMFD_Coupled {

 public:
  MatrixMFD_Coupled_TPFA(Teuchos::ParameterList& plist,
			 const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  MatrixMFD_Coupled_TPFA(const MatrixMFD_Coupled_TPFA& other);

  virtual void SetSubBlocks(const Teuchos::RCP<MatrixMFD>& blockA,
                    const Teuchos::RCP<MatrixMFD>& blockB);

  // Virtual copy constructor.
  virtual Teuchos::RCP<TreeMatrix> Clone() const {
    return Teuchos::rcp(new MatrixMFD_Coupled_TPFA(*this));
  }

  virtual int ApplyInverse(const TreeVector& X,
                            TreeVector& Y) const;
  virtual void ApplyInverse(const TreeVector& X,
                            const Teuchos::Ptr<TreeVector>& Y) const {
    ApplyInverse(X, *Y); }

  virtual void SymbolicAssembleGlobalMatrices();

 protected:

  virtual void AssembleSchur_() const;

  // sub-blocks
  Teuchos::RCP<MatrixMFD_TPFA> blockA_TPFA_;
  Teuchos::RCP<MatrixMFD_TPFA> blockB_TPFA_;

  Teuchos::RCP<Epetra_FEVbrMatrix> P2c2c_;

};

} // namespace
} // namespace

#endif
