#ifndef __AMANZI_CHEMISTRY_DIRECTSOLVER_HH__
#define __AMANZI_CHEMISTRY_DIRECTSOLVER_HH__

#include <iostream>
#include <vector>
using namespace std;

#include "chemistry/includes/solver.hh"
#include "chemistry/includes/Block.hpp"

class DirectSolver : public Solver {
  
public:
  DirectSolver(void);
  virtual ~DirectSolver();

  void Initialize(int n);
  void LUDecomposition(void);
  void LUBackSolve(void);
  void LUBackSolve(std::vector<double> &b);
  void Solve(void);
  void Solve(std::vector<double> &b);
  void Solve(MatrixBlock *A, std::vector<double> &b);
  void Solve(Block *A, std::vector<double> &b);
  
 private:

  double row_interchange_;
  std::vector<int> pivoting_indices_;
  std::vector<double> row_scaling_;
  bool factored_;

};

#endif // __AMANZI_CHEMISTRY_DIRECTSOLVER_HH__
