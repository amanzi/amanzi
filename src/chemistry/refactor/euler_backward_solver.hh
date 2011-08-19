/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_EULER_BACKWARD_SOLVER_HH_
#define AMANZI_CHEMISTRY_EULER_BACKWARD_SOLVER_HH_

//
// Class implementing a dedicated euler backward ode solver
//

#include <vector>
#include <ostream>

#include "ode_solver.hh"

namespace amanzi {
namespace chemistry {

class Block;

class EulerBackwardSolver : public ODESolver {
 public:
  EulerBackwardSolver();
  virtual ~EulerBackwardSolver();

  virtual void Solve(const std::vector<double>& func_eval, 
                     const Block& jacobian,
                     std::vector<double>* solution);
  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_EULER_BACKWARD_SOLVER_HH_
