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

  void Resize(int new_size);

  int size() const { return size_; };
  void set_size(int i) { size_ = i; };

  double* GetValues() const { return A_; };
  double& operator()(int i, int j) { return A_[j * size_ + i]; };
  const double& operator()(int i, int j) const { return A_[j * size_ + i]; };

  double GetRowAbsMax(int irow);

  // TODO(bandre): can we make some of these pointers const refs?
  void AddValue(int i, int j, double value);
  void AddValues(MatrixBlock* b);
  void AddValues(MatrixBlock* b, double scale);

  void ScaleRow(int irow, double scale);
  void ScaleColumn(int irow, double scale);
  void Scale(double scale);

  void Zero();
  void SetDiagonal(double d);

 private:
  int size_;
  double* A_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi

#endif
