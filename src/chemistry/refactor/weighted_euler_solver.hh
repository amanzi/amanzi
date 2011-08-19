/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_WEIGHTED_EULER_SOLVER_HH_
#define AMANZI_CHEMISTRY_WEIGHTED_EULER_SOLVER_HH_

//
// Class implementing the weighted euler ode solver, generalization of
// the crank-nicholson method that gives euler forward, euler backward
// and crank-nicholson depending on the value of a weighting
// parameter.
//
// df/dt = F(x,t)
//
// (f^{n+1} - f^{n})/dt = theta*F^{n+1} + (1.0 - theta)*F^{n}
// 
// when theta = 1.0, we get Euler backward
//      theta = 0.0, we get Euler Forward
//      theta = 0.5, we get Crank-Nicholson
//


#include <vector>
#include <ostream>

#include "ode_solver.hh"


namespace amanzi {
namespace chemistry {

class Block;

class WeightedEulerSolver : public ODESolver {
 public:
  WeightedEulerSolver(const double theta);
  virtual ~WeightedEulerSolver();

  virtual void Solve(const std::vector<double>& func_eval, 
                     const Block& jacobian,
                     std::vector<double>* solution);
  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_WEIGHTED_EULER_SOLVER_HH_
