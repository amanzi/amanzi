/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_ODE_SOLVER_HH_
#define AMANZI_CHEMISTRY_ODE_SOLVER_HH_

//
// Base class defining the API for ode solvers. possible inheriting
// classes for EulerForward, EulerBackward, RK, BDF, cvode wrappers,
// etc.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class Block;

class ODESolver {
 public:
  ODESolver();
  virtual ~ODESolver();

  virtual void Solve(const std::vector<double>& func_eval, 
                     const Block& jacobian,
                     std::vector<double>* solution) = 0;
  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_ODE_SOLVER_HH_
