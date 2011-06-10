/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_BLOCK_HH_
#define AMANZI_CHEMISTRY_BLOCK_HH_

// Boost may provide us with a more optimal matrix implementation - Glenn

class Block {
 public:
  Block();
  explicit Block(int n);
  virtual ~Block();

  int getSize(void) const {
    return this->size;
  };
  double** getValues(void) const {
    return this->A;
  };
  double GetValue(const int& i, const int& j) const {
    return this->A[i][j];
  };

  double getRowAbsMax(int imax);

  void setValue(int i, int j, double value);
  void setValues(double** values);
  void setValues(Block* b);
  void setValues(int ioffset, int joffset, Block* b);
  void setValues(double** values, double scale);
  void setValues(Block* b, double scale);
  void setValues(int ioffset, int joffset, Block* b, double scale);

  void addValue(int i, int j, double value);
  void addValues(double** values);
  void addValues(Block* b);
  void addValues(int ioffset, int joffset, Block* b);
  void addValues(double** values, double scale);
  void addValues(Block* b, double scale);
  void addValues(int ioffset, int joffset, Block* b, double scale);

  void scaleRow(int irow, double scale);
  void scaleColumn(int icol, double scale);
  void scale(double scale);

  void zero(void);
  void setDiagonal(double d);

  void print(void);


 private:

  int size;
  double** A;
};

#endif  // AMANZI_CHEMISTRY_BLOCK_HH_
