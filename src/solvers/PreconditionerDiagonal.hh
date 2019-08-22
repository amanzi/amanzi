/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PRECONDITIONER_DIAGONAL_HH_
#define AMANZI_PRECONDITIONER_DIAGONAL_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template<class Matrix, class Vector>
class PreconditionerDiagonal : public Preconditioner<Matrix,Vector> {
 public:
  PreconditionerDiagonal() {};
  ~PreconditionerDiagonal() {};

  void Init(const std::string& name, const Teuchos::ParameterList& plist) override {};
  void Update(const Teuchos::RCP<const Matrix>& A) override;
  void Destroy() override {};
  int ApplyInverse(const Vector& v, Vector& hv) const override;
  int returned_code() override { return 0; }

 private:
  Teuchos::RCP<Vector> work_vec_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi


//#ifdef HAVE_EPETRA_PRECONDITIONERS
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"

namespace Amanzi {
namespace AmanziPreconditioners {

template<>
void
PreconditionerDiagonal<Epetra_RowMatrix, Epetra_MultiVector>::Update(const Teuchos::RCP<const Epetra_RowMatrix>& matrix)
{
  work_vec_ = Teuchos::rcp(new Epetra_MultiVector(matrix->Map(), 1));
  matrix->ExtractDiagonalCopy(*(*work_vec_)(0));
}

template<>
int
PreconditionerDiagonal<Epetra_RowMatrix, Epetra_MultiVector>::ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) const
{
  return hv.ReciprocalMultiply(1.0, *work_vec_, v, 0.0); 
}

//#endif 
} // namespace
} // namespace

//#ifdef HAVE_TPETRA_PRECONDITIONERS

#include "VectorWrapper.hh"
#include "AmanziVector.hh"
#include "AmanziMatrix.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template<>
void
PreconditionerDiagonal<Matrix_type, VectorWrapper<Vector_type> >::Update(const Teuchos::RCP<const Matrix_type>& matrix)
{
  work_vec_ = Teuchos::rcp(new VectorWrapper<Vector_type>(Teuchos::rcp(new Vector_type(matrix->getRowMap()))));
  matrix->getLocalDiagCopy(*work_vec_->get());
}

template<>
int
PreconditionerDiagonal<Matrix_type, VectorWrapper<Vector_type> >::ApplyInverse(const VectorWrapper<Vector_type>& v, VectorWrapper<Vector_type>& hv) const
{
  return hv.ReciprocalMultiply(1.0, *work_vec_, v, 0.0); 
  return 0;
}

} // namespace
} // namespace
//#endif

#endif
