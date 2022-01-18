/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#ifndef AMANZI_CHEMISTRY_MATRIX_BLOCK_HH_
#define AMANZI_CHEMISTRY_MATRIX_BLOCK_HH_

#include "VerboseObject.hh"

#define PREFIX
#define F77_LAPACK_MANGLE(lcase,UCASE) lcase ## _
#define DGETRF_F77 F77_LAPACK_MANGLE(dgetrf,DGETRF)
#define DGETRI_F77 F77_LAPACK_MANGLE(dgetri,DGETRI)

#ifdef __cplusplus
extern "C" {
#endif

void PREFIX DGETRF_F77(int* nrow, int* ncol, double* a, int* lda, 
                       int* ipiv, int* info); 

void PREFIX DGETRI_F77(int* n, double* a, int* lda, 
                       int* ipiv, double* work, int* lwork, int* info); 

#ifdef __cplusplus
}
#endif


namespace Amanzi {
namespace AmanziChemistry {

class MatrixBlock {
 public:
  MatrixBlock();
  explicit MatrixBlock(int n);
  MatrixBlock(int n, int m);
  ~MatrixBlock();

  void Resize(int new_size);

  int size() const { return size_; };
  int cols() const { return cols_; };

  double* GetValues() const { return A_; };
  double& operator()(int i, int j) { return A_[j * size_ + i]; };
  const double& operator()(int i, int j) const { return A_[j * size_ + i]; };

  double GetRowAbsMax(int irow);

  void AddValue(int i, int j, double value);
  void AddValues(MatrixBlock* b);
  void AddValues(MatrixBlock* b, double scale);

  void ScaleRow(int irow, double scale);
  void ScaleColumn(int irow, double scale);
  void Scale(double scale);

  void Zero();
  void SetDiagonal(double d);

  // non trivilar routines
  int Inverse();

  // output 
  friend std::ostream& operator << (std::ostream& os, const MatrixBlock& A) {
    for (int i = 0; i < A.size(); i++) {
      for (int j = 0; j < A.cols(); j++) {
        os << std::setw(12) << std::setprecision(12) << A(i, j) << " ";
      }
      os << "\n";
    }
    return os;
  }

 private:
  int size_, cols_;
  double* A_;
};

// non-member functions
void Multiply(const MatrixBlock& A, const MatrixBlock& B, MatrixBlock& AB);

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
