/*
This is the Audit component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef __MATRIX_AUDIT_HPP__
#define __MATRIX_AUDIT_HPP__

#include "Epetra_Operator.h"

#include "Matrix_MFD.hpp"

namespace Amanzi {

int const MATRIX_AUDIT_MFD = 1;

class Matrix_Audit {
 public:
  Matrix_Audit(Teuchos::RCP<AmanziMesh::Mesh> mesh, AmanziFlow::Matrix_MFD* matrix);
  ~Matrix_Audit();

  // main members
  void InitAudit();
  int CheckSpectralBounds();

 private:
  void OrderByIncrease(int n, double* mem);

  int matrix_type;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  AmanziFlow::Matrix_MFD* matrix_;

  std::vector<Teuchos::SerialDenseMatrix<int, double> >* A;  // local matrices
  int lda;  // maximal size of elemental matrices

  int lwork1;  // work memory
  double *dmem1, *dmem2;
  double *dwork1;
};

}  // namespace Amanzi

#endif
