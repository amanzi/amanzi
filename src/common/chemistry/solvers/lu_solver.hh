/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Simple class template
*/

#ifndef AMANZI_CHEMISTRY_LU_SOLVER_HH_
#define AMANZI_CHEMISTRY_LU_SOLVER_HH_

#include <ostream>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

class MatrixBlock;

class LUSolver {
 public:
  LUSolver();
  virtual ~LUSolver() {};

  void Initialize(const int size);
  void Solve(MatrixBlock* A, std::vector<double>* b);

  static const double kSmallNumber;

 protected:
  void Decomposition(MatrixBlock* A);
  void BackSolve(MatrixBlock* A, std::vector<double>* b);
  void set_system_size(const int size) {
    this->system_size_ = size;
  }
  int system_size(void) const {
    return this->system_size_;
  }

 private:
  int system_size_;
  double row_interchange_;
  std::vector<int> pivoting_indices_;
  std::vector<double> row_scaling_;
  bool factored_;
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
