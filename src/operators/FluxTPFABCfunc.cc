#include "OperatorDefs.hh"
#include "OperatorDiffusionFV.hh"
#include "FluxTPFABCfunc.hh"

namespace Amanzi {
namespace Operators {


template <class WRM>
double solve_nonlinear_bisection(double face_val, 
                                 double atm_press, 
                                 FlowFluxTPFA_BCFunc<WRM>& func,
                                 Tol_& tol,
                                 int max_it,
                                 int& actual_it){

  double res = func(face_val);
  double left = 0.;
  double right = 0.;
  double lres = 0.;
  double rres = 0.;

  if (res > 0.) {
    left = face_val;
    lres = res;
    right = std::max(face_val, atm_press);
    rres = func(right);
    while (rres > 0.) {
      right += atm_press;
      rres = func(right);
    }
  }
  else {
    right = face_val;
    rres = res;
#if DEBUG_FLAG
    std::cout << "RIGHT = " << right << ", " << rres << std::endl;
#endif
    left = std::min(atm_press, face_val);
    lres = func(left);
    while (lres < 0.) {
#if DEBUG_FLAG
      std::cout << "LEFT = " << left << ", " << lres << std::endl;
#endif
      left -= atm_press;
      lres = func(left);
    }
  }
#if DEBUG_FLAG
  std::cout << "   bracket (res): " << left << " (" << lres << "), "
            << right << " (" << rres << ")" << std::endl;
#endif

  std::pair<double,double> result =
    boost::math::tools::toms748_solve(func, left, right, lres, rres, tol, actual_it);
  if (actual_it >= max_it) {
    std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
    return 3;
  }
	  
  face_val = (result.first + result.second) / 2.;

#if DEBUG_FLAG
  std::cout << "face_val = "<<face_val<<"\n";
#endif
  return face_val;


}


}
}
