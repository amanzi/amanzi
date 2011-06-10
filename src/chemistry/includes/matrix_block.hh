/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_MATRIXBLOCK_HH_
#define AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_MATRIXBLOCK_HH_

#include "chemistry/includes/block.hh"

// Boost may provide us with a more optimal matrix implementation - Glenn

class MatrixBlock {

 public:
  MatrixBlock();
  MatrixBlock(int);
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

  double GetRowAbsMax(int);

  void SetValue(int, int, double);
  void SetValues(double**);
  void SetValues(MatrixBlock*);
  void SetValues(Block*);
  void SetValues(int, int, MatrixBlock*);
  void SetValues(double**, double scale);
  void SetValues(MatrixBlock*, double scale);
  void SetValues(int, int, MatrixBlock*, double scale);

  void AddValue(int, int, double);
  void AddValues(double**);
  void AddValues(MatrixBlock*);
  void AddValues(int, int, MatrixBlock*);
  void AddValues(double**, double scale);
  void AddValues(MatrixBlock*, double scale);
  void AddValues(int, int, MatrixBlock*, double scale);

  void ScaleRow(int, double);
  void ScaleColumn(int, double);
  void Scale(double);

  void Zero(void);
  void SetDiagonal(double d);

  void Print(void);


 private:

  int size_;
  double** A_;

};

#endif  // AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_MATRIXBLOCK_HH_
