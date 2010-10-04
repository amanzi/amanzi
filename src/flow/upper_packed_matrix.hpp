#ifndef __UPPER_PACKED_MATRIX_HPP__
#define __UPPER_PACKED_MATRIX_HPP__

class upper_packed_matrix {

public:

  upper_packed_matrix(double* a_, int n_) : a(a_), n(n_) {};
  ~upper_packed_matrix() {};

  void invert();
  void factor();
  void solve(double * x);
  void col_sum(double * csum);
  void sym_matmul(double* xin, double* xresult);

private:
  double* a;
  int n;


};

#endif


