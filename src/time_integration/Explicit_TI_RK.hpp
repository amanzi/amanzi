#ifndef _EXPLICIT_TI_RK_HPP_
#define _EXPLICIT_TI_RK_HPP_

#include "Epetra_Vector.h"

#include "Explicit_TI_fnBase.hpp"
#include "boost/numeric/ublas/matrix.hpp"

namespace Amanzi {
namespace Explicit_TI {

class RK {
  // this class implements several explicit Runge Kutta methods:
  // forward Euler     (1st order)  --> forward_euler
  // Heun-Euler method (2nd order)  --> heun_euler 
  // Midpoint method   (2nd order)  --> midpoint 
  // Ralston method    (2nd order)  --> ralston 
  // Kutta method      (3rd order)  --> kutta_3rd_order
  // Runge Kutta       (4th order)  --> runge_kutta_4th_order
  // User defined      (whatever)   --> user_defined, use special constructor to create

  // the RK tableau is made up of the three private objects a, b, and c below.
  // they are arranged as follows:
  //
  // c[0]   | 
  // c[1]   | a(1,0)
  //  .     | a(2,0)    a(2,1)
  //  .     |   .         .        
  //  .     |   .         . 
  // c[s-1 ]| a(s-1,0)  a(s-1,1)  . . .  a(s-1,s-2) 
  // ---------------------------------------------------------
  //        |   b[0]      b[1]    . . .    b[s-2]      b[s-1] 
  // 
  // note that c[0] should always equal zero, and that the entries in the matrix
  // a that are not listed in this tableau are not used
  //
  // the implemented general Runge Kutta scheme of order s based on this tableau arrangement is
  //
  //     y_{n+1} = y_n + \sum{i=0}^{s-1} b[i]*k_i
  //
  // with 
  //
  //     k_0 = h * f(t_n, y_n)
  //     k_1 = h * f(t_n + c[1]*h, y_n + a(1,0)*k_0)
  //     k_2 = h * f(t_n + c[2]*h, y_n + a(2,0)*k_0 + a(2,1)*k_1)
  //      .
  //      .
  //      .
  //     k_{s-1} = h * f(t_n + c[s-1]*h, y_n + a(s-1,0)*k_0 + ... + a(s-1,s-2)*k_{s-2})
 
 public:
  enum method_t{forward_euler, 
	        heun_euler, 
                midpoint, 
                ralston, 
                kutta_3rd_order, 
                runge_kutta_4th_order,
                user_defined};

  // constructor for pre-coded RK methods (see list of methods in the method_t type above)
  RK(Explicit_TI::fnBase& fn_, 
     const method_t method, 
     const Epetra_Vector& example_vector_);

  // constructor for user defined RK methods
  RK(Explicit_TI::fnBase& fn_,
     const int order_,
     const boost::numeric::ublas::matrix<double> a_,
     const std::vector<double> b_,
     const std::vector<double> c_,
     const Epetra_Vector& example_vector_);
      
  ~RK();

  void step(const double t, const double h, const Epetra_Vector& y, Epetra_Vector& y_new);

  int get_order() {return order;}
  method_t get_method() {return method;}

 private:
  void set_method(const Explicit_TI::RK::method_t method_);
  void create_storage(const Epetra_Vector& example_vector_); 
  void delete_storage();

  Explicit_TI::fnBase& fn;
  int order;
  method_t method;
  boost::numeric::ublas::matrix<double> a;
  std::vector<double> b, c;

  std::vector<Epetra_Vector*> k;
};

}  // namespace Explicit_TI
}  // namespace Amanzi


#endif // _EXPLICIT_TI_RK_HPP_
