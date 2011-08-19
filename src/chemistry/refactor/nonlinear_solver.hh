/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_NONLINEAR_SOLVER_HH_
#define AMANZI_CHEMISTRY_NONLINEAR_SOLVER_HH_

//
// Base class defining the API for all nonlinear solver classes.
// Inheriting classes would be newton-raphson, newton-raphson with
// under relaxation, newton-raphson with line search, etc. Could also
// have a wrapper class for kinsol.
//

#include <vector>
#include <ostream>


namespace amanzi {
namespace chemistry {

class NonlinearSolver {
 public:
  NonlinearSolver();
  virtual ~NonlinearSolver();

  virtual void Display(std::ostream& output) const;

 protected:
 private:

};

}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_NONLINEAR_SOLVER_HH_
