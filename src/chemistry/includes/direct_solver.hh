/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_DIRECTSOLVER_HH_
#define AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_DIRECTSOLVER_HH_

#include <iostream>
#include <vector>

#include "chemistry/includes/solver.hh"
#include "chemistry/includes/block.hh"

class DirectSolver : public Solver {
 public:
  DirectSolver(void);
  virtual ~DirectSolver();

  void Initialize(int n);
  void LUDecomposition(void);
  void LUBackSolve(void);
  void LUBackSolve(std::vector<double>* b);
  void Solve(void);
  void Solve(std::vector<double>* b);
  void Solve(MatrixBlock* A, std::vector<double>* b);
  void Solve(Block* A, std::vector<double>* b);

 private:

  double row_interchange_;
  std::vector<int> pivoting_indices_;
  std::vector<double> row_scaling_;
  bool factored_;
};

#endif  // AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_DIRECTSOLVER_HH_
