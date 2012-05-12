/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_LU_SOLVER_HH_
#define AMANZI_CHEMISTRY_LU_SOLVER_HH_

//
// Simple class template
//

#include <vector>
#include <ostream>

#include "linear_solver.hh"

namespace amanzi {
namespace chemistry {

class MatrixBlock;

class LUSolver : public LinearSolver {
 public:
  LUSolver();
  virtual ~LUSolver();

  void Initialize(const int size);
  void Solve(MatrixBlock* A, std::vector<double>* b);

  static const double kSmallNumber;

 protected:
  void Decomposition(MatrixBlock* A);
  void BackSolve(MatrixBlock* A, std::vector<double>* b);
 private:
  double row_interchange_;
  std::vector<int> pivoting_indices_;
  std::vector<double> row_scaling_;
  bool factored_;
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_LU_SOLVER_HH_
