/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_NEWTON_RAPHSON_SOLVER_HH_
#define AMANZI_CHEMISTRY_NEWTON_RAPHSON_SOLVER_HH_

//
// Implement Newton-Raphson algorithm (with under relaxation?) for
// solution of a system of linear equations. NR with line search, etc
// would be a different subclass with the same interface.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class NewtonRaphsonSolver : public NonlinearSolver {
 public:
  NewtonRaphsonSolver();
  virtual ~NewtonRaphsonSolver();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_NEWTON_RAPHSON_SOLVER_HH_
