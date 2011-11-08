#include "UnitTest++.h"

#include "Explicit_TI_fnBase.hpp"
#include "Explicit_TI_RK.hpp"

#include "Epetra_BlockMap.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialComm.h"

SUITE(TimeIntegrationTests) {
using namespace Amanzi;

  // ODE: y' = y
  class fn1 : public Explicit_TI::fnBase {
  public:
    void fun(const double t, const Epetra_Vector& y, Epetra_Vector& y_new)
    {
      y_new = y;
    }
  };


  TEST(Explicit_RK_Euler) {
    cout << "Test: Explicit_RK_Euler" << endl;    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    Explicit_TI::RK::method_t method = Explicit_TI::RK::forward_euler;
    Explicit_TI::RK explicit_time_integrator(f, method, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t,h,y,y_new);
	t=t+h;
	y = y_new;
      }
    while (t<1.0);
    CHECK_CLOSE(y[0],exp(t),2.0*h);
  }
       

  TEST(Explicit_RK_Heun) {
    cout << "Test: Explicit_RK_Heun" << endl;    
    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    Explicit_TI::RK::method_t method = Explicit_TI::RK::heun_euler;
    Explicit_TI::RK explicit_time_integrator(f, method, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t,h,y,y_new);
	t=t+h;
	y = y_new;
      }
    while (t<1.0);
    CHECK_CLOSE(y[0],exp(t),pow(h,2));
  }
       

  TEST(Explicit_RK_Midpoint) {
    cout << "Test: Explicit_RK_Midpoint" << endl;    
    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    Explicit_TI::RK::method_t method = Explicit_TI::RK::midpoint;
    Explicit_TI::RK explicit_time_integrator(f, method, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t,h,y,y_new);
	t=t+h;
	y = y_new;
      }
    while (t<1.0);
    CHECK_CLOSE(y[0],exp(t),pow(h,2));
  }

  TEST(Explicit_RK_Ralston) {
    cout << "Test: Explicit_RK_Rapson" << endl;    
    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    Explicit_TI::RK::method_t method = Explicit_TI::RK::ralston;
    Explicit_TI::RK explicit_time_integrator(f, method, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t,h,y,y_new);
	t=t+h;
	y = y_new;
      }
    while (t<1.0);
    CHECK_CLOSE(y[0],exp(t),pow(h,2));
  }


  TEST(Explicit_RK_Kutta3D) {
    cout << "Test: Explicit_RK_Kutta3D" << endl;    
    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    Explicit_TI::RK::method_t method = Explicit_TI::RK::kutta_3rd_order;
    Explicit_TI::RK explicit_time_integrator(f, method, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t,h,y,y_new);
	t=t+h;
	y = y_new;
      }
    while (t<1.0);
    CHECK_CLOSE(y[0],exp(t),pow(h,3));
  }

  TEST(Explicit_RK_UserDefined) {
    cout << "Test: Explicit_RK_UserDefined" << endl;    
    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    int order = 2;
    boost::numeric::ublas::matrix<double> a(order,order);
    std::vector<double> b(order);
    std::vector<double> c(order);
    
    a(1,0) = 1.0;

    b[0] = 0.5;
    b[1] = 0.5;

    c[0] = 0.0;
    c[1] = 1.0;

    Explicit_TI::RK explicit_time_integrator(f, order, a, b, c, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t,h,y,y_new);
	t = t + h;
	y = y_new;
      }
    while (t<1.0);


    CHECK_CLOSE(y[0],exp(t),pow(h,2));
  }

  TEST(Explicit_RK_RK4) {
    cout << "Test: Explicit_RK_RK4" << endl;    
    
    Epetra_Comm* comm = new Epetra_SerialComm();    
    Epetra_BlockMap map(1,1,0,*comm);
    Epetra_Vector y(map);
    Epetra_Vector y_new(map);

    fn1 f;
    Explicit_TI::RK::method_t method = Explicit_TI::RK::runge_kutta_4th_order;
    Explicit_TI::RK explicit_time_integrator(f, method, y); 
		
    // initial value
    y.PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;
    
    // integrate to t=1.0
    do 
      {
	explicit_time_integrator.step(t, h, y, y_new);
	t = t + h;
	y = y_new;
      }
    while (t<1.0);
    CHECK_CLOSE(y[0],exp(t),pow(h,4));
  }
}
