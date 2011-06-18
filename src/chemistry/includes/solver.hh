/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_SOLVER_HH_
#define AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_SOLVER_HH_

#include <vector>

#include "chemistry/includes/matrix_block.hh"

namespace amanzi {
namespace chemistry {

class Solver {
 public:
  Solver();
  virtual ~Solver();

  virtual void Initialize(int i);
  virtual void Solve(void) = 0;

  void set_system_size(int i) {
    this->system_size_ = i;
  }
  int system_size(void) const {
    return this->system_size_;
  }

 protected:

  int system_size_;
  std::vector<double> right_hand_side_;
  std::vector<double> solution_;
  MatrixBlock* A_;
};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_AMANZI_CHEMISTRY_SOLVER_HH_
