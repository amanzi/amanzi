#include "ExplicitTIRK.hh"

namespace Amanzi {

ExplicitTIRK::ExplicitTIRK(ExplicitTIfnBase& fn_, 
       const ExplicitTIRK::method_t method_, 
       const TreeVector& example_vector):
    fn(fn_)
{ 
  set_method(method_);
  create_storage(example_vector); 
}


ExplicitTIRK::ExplicitTIRK(ExplicitTIfnBase& fn_,
       const int order_,
       const boost::numeric::ublas::matrix<double> a_,
       const std::vector<double> b_,
       const std::vector<double> c_,
       const TreeVector& example_vector):
  fn(fn_), order(order_), method(ExplicitTIRK::user_defined)
{
  a.resize(a_.size1(), a_.size2());
  b.resize(b_.size());
  c.resize(c_.size());

  a = a_;
  b = b_;
  c = c_;

  create_storage(example_vector); 
}


ExplicitTIRK::~ExplicitTIRK()
{
  delete_storage();
}


void ExplicitTIRK::set_method(const ExplicitTIRK::method_t method_)
{
  method = method_;
  
  switch (method) {
  case forward_euler: 
    order = 1;
    break;
  case heun_euler:
    order = 2;
    break;
  case midpoint:
    order = 2;
    break;
  case ralston:
    order = 2;
    break;
  case kutta_3rd_order:
    order = 3;
    break;
  case runge_kutta_4th_order:
    order = 4;
    break;
  }

  a.resize(order, order);
  b.resize(order);
  c.resize(order);

  // initialize the RK tableau
  switch (method) {

  case forward_euler:
    b[0] = 1.0;
    c[0] = 0.0;
    break;

  case heun_euler:
    a(1,0) = 1.0;

    b[0] = 0.5;
    b[1] = 0.5;

    c[0] = 0.0;
    c[1] = 1.0;
    break;

  case midpoint:
    a(1,0) = 0.5;
    
    b[0] = 0.0;
    b[1] = 1.0;

    c[0] = 0.0;
    c[1] = 0.5;
    break;

  case ralston:
    a(1,0) = 2.0/3.0;
    
    b[0] = 0.25;
    b[1] = 0.75;
    
    c[0] = 0.0;
    c[1] = 2.0/3.0;
    break;

  case kutta_3rd_order:
    a(1,0) = 0.5;
    a(2,0) = -1.0;
    
    a(2,1) = 2.0;

    b[0] = 1.0/6.0;
    b[1] = 2.0/3.0;
    b[2] = 1.0/6.0;

    c[0] = 0.0;
    c[1] = 0.5;
    c[2] = 1.0;
    break;

  case runge_kutta_4th_order:
    a(1,0) = 0.5;
    a(2,0) = 0.0;
    a(3,0) = 0.0;
    
    a(2,1) = 0.5;
    a(3,1) = 0.0;
    
    a(3,2) = 1.0;
    
    b[0] = 1.0/6.0;
    b[1] = 1.0/3.0;
    b[2] = 1.0/3.0;
    b[3] = 1.0/6.0;

    c[0] = 0.0;
    c[1] = 0.5;
    c[2] = 0.5;
    c[3] = 1.0;
  }
}


void ExplicitTIRK::create_storage(const TreeVector& example_vector)
{
  k.resize(order);
  for (int i=0; i<order; i++)
    {
      k[i] = new TreeVector(example_vector);
    }
}


void ExplicitTIRK::delete_storage()
{
  for (int i=0; i<k.size(); i++)
    {
      delete k[i];
    }
}


void ExplicitTIRK::step(const double t, const double h, const TreeVector& y, TreeVector& y_new)
{
  TreeVector sum_vec(y);
  double sum_time;
  for (int i=0; i<order; i++) 
    {
      sum_time = t + c[i]*h;
      
      sum_vec = y;
      
      for (int j=0; j<i; j++) 
        {
          if (a(i,j) != 0.0) 
	    {
	      sum_vec.Update(a(i,j), *k[j], 1.0);
	    }
	}
      fn.fun(sum_time, sum_vec, *k[i]);
      k[i]->Scale(h);
    }

  y_new = y;
      
  for (int i=0; i<order; i++)
    {
      if (b[i] != 0.0) 
        {
           y_new.Update(b[i],*k[i],1.0);
        }
    }
}

}  // namespace Amanzi
