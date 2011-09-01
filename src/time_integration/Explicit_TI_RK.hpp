#ifndef _EXPLICIT_TI_RK_HPP_
#define _EXPLICIT_TI_RK_HPP_

#include "Epetra_MultiVector.h"

#include "Explicit_TI_fnBase.hpp"
#include "boost/numeric/ublas/matrix.hpp"

namespace Explicit_TI {



  class RK {

    // this class implements several explicit Runge Kutta methods:
    // forward Euler (1st order)
    // Heun's method (2nd order)
    // Kutta method  (3rd order)
    // Runge Kutta   (4th order)

  public:

    enum method_t { forward_euler, heun_method, kutta_3rd_order, runge_kutta_4th_order };

    RK(Explicit_TI::fnBase& fn_, const method_t method, const Epetra_MultiVector& example_vector_);
    ~RK();

    void step(const double t, const double h, const Epetra_MultiVector& y, Epetra_MultiVector& y_new);

    int get_order () { return order; }
    method_t get_method () { return method; }

  private:
    void set_method ( const Explicit_TI::RK::method_t method_ );
    void create_storage ( const Epetra_MultiVector& example_vector_ ); 
    void delete_storage ( );

  private:
    Explicit_TI::fnBase& fn;
    int order;
    method_t method;
    boost::numeric::ublas::matrix<double> a;
    std::vector<double> b, c;

    std::vector<Epetra_MultiVector*> k;
  };

} // namespace Explicit_TI





#endif // _EXPLICIT_TI_RK_HPP_
