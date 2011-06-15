/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_MATRIXBLOCK_HH_
#define AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_MATRIXBLOCK_HH_

#include "chemistry/includes/block.hh"

// Boost may provide us with a more optimal matrix implementation - Glenn

class MatrixBlock {
 public:
  MatrixBlock();
  explicit MatrixBlock(int n);
  virtual ~MatrixBlock();

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

  void SetValue(int i, int j, double value);
  void SetValues(double** values);
  void SetValues(MatrixBlock* b);
  void SetValues(Block* b);
  void SetValues(int ioffset, int ioffset, MatrixBlock* b);
  void SetValues(double** values, double scale);
  void SetValues(MatrixBlock* b, double scale);
  void SetValues(int ioffset, int ioffset, MatrixBlock* b, double scale);

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

  void Print(void);


 private:

  int size_;
  double** A_;
};

#endif  // AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_MATRIXBLOCK_HH_
