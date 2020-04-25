/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Markus Berndt (berndt@lanl.gov)
*/
//! Explicit time integration methods in a generalized form.

/*!

This class implements several explicit Runge Kutta methods:

- forward Euler     (1st order)  --> "forward_euler"
- Heun-Euler method (2nd order)  --> "heun_euler"
- Midpoint method   (2nd order)  --> "midpoint"
- Ralston method    (2nd order)  --> "ralston"
- TVD RK method     (3rd order)  --> "tvd_3rd_order"
- Kutta method      (3rd order)  --> "kutta_3rd_order"
- Runge Kutta       (4th order)  --> "runge_kutta_4th_order"
- User defined      (whatever)   --> user_defined, use special constructor to create

Note that user-defined is only for developers currently, and cannot be created from an input file.

The RK tableau is made up of the three private objects a, b, and c below.  they
are arranged as follows:

.. code-block:

    c[0]   | 
    c[1]   | a(1,0)
     .     | a(2,0)    a(2,1)
     .     |   .         .        
     .     |   .         . 
    c[s-1 ]| a(s-1,0)  a(s-1,1)  . . .  a(s-1,s-2) 
    ---------------------------------------------------------
           |   b[0]      b[1]    . . .    b[s-2]      b[s-1] 
  
Note that c[0] should always equal zero, and that the entries in the matrix a
that are not listed in this tableau are not used
  
The implemented general Runge Kutta scheme of order s based on this tableau arrangement is

.. math::
    y_{n+1} = y_n + \sum{i=0}^{s-1} b[i]*k_i
  
    with 
  
      k_0 = h * f(t_n, y_n) \\
      k_1 = h * f(t_n + c[1]*h, y_n + a(1,0)*k_0) \\
      k_2 = h * f(t_n + c[2]*h, y_n + a(2,0)*k_0 + a(2,1)*k_1) \\
       . \\
       . \\
       . \\ 
      k_{s-1} = h * f(t_n + c[s-1]*h, y_n + a(s-1,0)*k_0 + ... + a(s-1,s-2)*k_{s-2})


.. _explicit-ti-rk-spec:
.. admonition:: explicit-ti-rk-spec

    * `"verbose object`" ``[verbose-object-spec]`` A `Verbose Object`_

    * `"RK method`" ``[string]`` **forward euler**  One of: `"forward Euler`", `"heun euler`", `"midpoint`", `"ralston`", `"tvd 3rd order`", `"kutta 3rd order`", `"runge kutta 4th order`"
      
*/



#ifndef AMANZI_EXPLICIT_TI_RK_HH_
#define AMANZI_EXPLICIT_TI_RK_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"

#include "errors.hh"
#include "Explicit_TI_FnBase.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Explicit_TI {

enum method_t{forward_euler, 
              heun_euler, 
              midpoint, 
              ralston, 
              tvd_3rd_order,
              kutta_3rd_order, 
              runge_kutta_4th_order,
              user_defined};

template<class Vector>
class RK {
 
 public:
  // constructor for pre-coded RK methods (see list of methods in the method_t type above)
  RK(fnBase<Vector>& fn, 
     const method_t method, 
     const Vector& initvector);

  // constructor from ParameterList
  RK(fnBase<Vector>& fn,
     Teuchos::ParameterList& plist,
     const Vector& initvector);
     
  // constructor for user defined RK methods
  RK(fnBase<Vector>& fn,
     const int order,
     const Epetra_SerialDenseMatrix& a,
     const std::vector<double> b,
     const std::vector<double> c,
     const Vector& initvector);

  void TimeStep(const double t, const double h, const Vector& y, Vector& y_new);

  int order() { return order_; }
  method_t method() { return method_; }

 private:
  void InitMethod_(const method_t method);
  void CreateStorage_(const Vector& initvector);

  fnBase<Vector>& fn_;
  int order_;
  method_t method_;
  Epetra_SerialDenseMatrix a_;
  std::vector<double> b_, c_;
  std::vector<Teuchos::RCP<Vector> > k_;

  Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;
};


template<class Vector>
RK<Vector>::RK(fnBase<Vector>& fn, 
               const method_t method, 
               const Vector& initvector) :
    fn_(fn)
{ 
  InitMethod_(method);
  CreateStorage_(initvector);
  vo_ = Teuchos::rcp(new VerboseObject("TI::RK", plist_));
}


template<class Vector>
RK<Vector>::RK(fnBase<Vector>& fn, 
               Teuchos::ParameterList& plist,
               const Vector& initvector) :
    fn_(fn),
    plist_(plist)
{ 
  std::string methodstring = plist_.get<std::string>("RK method", "forward euler");
  method_t method;
  if (methodstring == "forward euler") {
    method = forward_euler;
  } else if (methodstring == "heun euler") {
    method = heun_euler;
  } else if (methodstring == "midpoint") {
    method = midpoint;
  } else if (methodstring == "ralston") {
    method = ralston;
  } else if (methodstring == "tvd 3rd order") {
    method = tvd_3rd_order;
  } else if (methodstring == "kutta 3rd order") {
    method = kutta_3rd_order;
  } else if (methodstring == "runge kutta 4th order") {
    method = runge_kutta_4th_order;
  } else {
    Errors::Message msg("RK TI: wrong value of parameter `\"RK method`\"");
    Exceptions::amanzi_throw(msg);
  }

  // should add way to specify user defined method from PList. --etc

  InitMethod_(method);
  CreateStorage_(initvector); 
  vo_ = Teuchos::rcp(new VerboseObject("TI::RK",plist_));
}


template<class Vector>
RK<Vector>::RK(fnBase<Vector>& fn,
               const int order,
               const Epetra_SerialDenseMatrix& a,
               const std::vector<double> b,
               const std::vector<double> c,
               const Vector& initvector) :
    fn_(fn),
    order_(order),
    method_(user_defined),
    a_(a)
{
  b_.resize(b.size());
  c_.resize(c.size());
  b_ = b;
  c_ = c;

  CreateStorage_(initvector); 
  vo_ = Teuchos::rcp(new VerboseObject("TI::RK",plist_));
}


template<class Vector>
void RK<Vector>::InitMethod_(const method_t method)
{
  method_ = method;
  
  switch (method_) {
  case forward_euler: 
    order_ = 1;
    break;
  case heun_euler:
    order_ = 2;
    break;
  case midpoint:
    order_ = 2;
    break;
  case ralston:
    order_ = 2;
    break;
  case tvd_3rd_order:
    order_ = 3;
    break;
  case kutta_3rd_order:
    order_ = 3;
    break;
  case runge_kutta_4th_order:
    order_ = 4;
    break;
  default:
    order_ = -1;
  }

  a_.Shape(order_, order_);
  b_.resize(order_);
  c_.resize(order_);

  // initialize the RK tableau
  switch (method) {

  case forward_euler:
    b_[0] = 1.0;
    c_[0] = 0.0;
    break;

  case heun_euler:
    a_(1,0) = 1.0;

    b_[0] = 0.5;
    b_[1] = 0.5;

    c_[0] = 0.0;
    c_[1] = 1.0;
    break;

  case midpoint:
    a_(1,0) = 0.5;
    
    b_[0] = 0.0;
    b_[1] = 1.0;

    c_[0] = 0.0;
    c_[1] = 0.5;
    break;

  case ralston:
    a_(1,0) = 2.0/3.0;
    
    b_[0] = 0.25;
    b_[1] = 0.75;
    
    c_[0] = 0.0;
    c_[1] = 2.0/3.0;
    break;

  case tvd_3rd_order:
    a_(1,0) = 1.0;
    a_(2,0) = 0.25;
    
    a_(2,1) = 0.25;

    b_[0] = 1.0/6.0;
    b_[1] = 1.0/6.0;
    b_[2] = 2.0/3.0;

    c_[0] = 0.0;
    c_[1] = 1.0;
    c_[2] = 0.5;
    break;

  case kutta_3rd_order:
    a_(1,0) = 0.5;
    a_(2,0) = -1.0;
    
    a_(2,1) = 2.0;

    b_[0] = 1.0/6.0;
    b_[1] = 2.0/3.0;
    b_[2] = 1.0/6.0;

    c_[0] = 0.0;
    c_[1] = 0.5;
    c_[2] = 1.0;
    break;

  case runge_kutta_4th_order:
    a_(1,0) = 0.5;
    a_(2,0) = 0.0;
    a_(3,0) = 0.0;
    
    a_(2,1) = 0.5;
    a_(3,1) = 0.0;
    
    a_(3,2) = 1.0;
    
    b_[0] = 1.0/6.0;
    b_[1] = 1.0/3.0;
    b_[2] = 1.0/3.0;
    b_[3] = 1.0/6.0;

    c_[0] = 0.0;
    c_[1] = 0.5;
    c_[2] = 0.5;
    c_[3] = 1.0;
    break;
    default:
      break;
  }
}


template<class Vector>
void RK<Vector>::CreateStorage_(const Vector& initvector)
{
  k_.resize(order_);
  for (int i=0; i!=order_; ++i) {
    k_[i] = Teuchos::rcp(new Vector(initvector));
  }
}


template<class Vector>
void RK<Vector>::TimeStep(double t, double h, const Vector& y, Vector& y_new)
{
  Vector y_tmp(y);
  fn_.ModifySolution(t, y_tmp);

  double sum_time;
  for (int i = 0; i != order_; ++i) {
    sum_time = t + c_[i] * h;

    if (i == 0) {
      fn_.FunctionalTimeDerivative(sum_time, y_tmp, *k_[0]);
    } else {
      y_new = y_tmp;
      
      for (int j = 0; j != i; ++j) {
        if (a_(i,j) != 0.0) {
          y_new.Update(a_(i,j), *k_[j], 1.0);
        }
      }
      fn_.ModifySolution(sum_time, y_new);
      fn_.FunctionalTimeDerivative(sum_time, y_new, *k_[i]);
    }

    k_[i]->Scale(h);
  }

  y_new = y_tmp;
  for (int i = 0; i != order_; ++i) {
    if (b_[i] != 0.0) {
      y_new.Update(b_[i], *k_[i], 1.0);
    }
  }
}


}  // namespace Explicit_TI
}  // namespace Amanzi


#endif // AMANZI_EXPLICIT_TI_RK_HPP_
