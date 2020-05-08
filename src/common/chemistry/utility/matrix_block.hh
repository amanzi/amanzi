/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre
*/

#ifndef AMANZI_CHEMISTRY_MATRIXBLOCK_HH_
#define AMANZI_CHEMISTRY_MATRIXBLOCK_HH_


// Boost may provide us with a more optimal matrix implementation - Glenn

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

class MatrixBlock {
 public:
  MatrixBlock();
  explicit MatrixBlock(int n);
  virtual ~MatrixBlock();

  void Resize(const int new_size);

  int size(void) const {
    return this->size_;
  };
  void set_size(int i) {
    this->size_ = i;
  };

  double** GetValues(void) const {
    return this->A_;
  };
  double** GetValuesMutable(void) {
    return this->A_;
  };
  double GetValue(const int& i, const int& j) const {
    return this->A_[i][j];
  };

  double GetRowAbsMax(int irow);

  // TODO(bandre): can we make some of these pointers const refs?
  void SetValue(int i, int j, double value);
  void SetValues(double** values);
  void SetValues(MatrixBlock* b);
  void SetValues(int ioffset, int joffset, MatrixBlock* b);
  void SetValues(double** values, double scale);
  void SetValues(MatrixBlock* b, double scale);
  void SetValues(int ioffset, int joffset, MatrixBlock* b, double scale);

  void AddValue(int i, int j, double value);
  void AddValues(double** values);
  void AddValues(MatrixBlock* b);
  void AddValues(int ioffset, int joffset, MatrixBlock* b);
  void AddValues(double** values, double scale);
  void AddValues(MatrixBlock* b, double scale);
  void AddValues(int ioffset, int joffset, MatrixBlock* b, double scale);

  void ScaleRow(int irow, double scale);
  void ScaleColumn(int irow, double scale);
  void Scale(double scale);

  void Zero(void);
  void SetDiagonal(double d);

 private:
  void AllocateMemory(void);
  void FreeMemory(void);

  int size_;
  double** A_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
