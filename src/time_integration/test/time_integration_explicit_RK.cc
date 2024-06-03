/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include "UnitTest++.h"

#include "Teuchos_SerialDenseMatrix.h"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"

#include "Explicit_TI_FnBase.hh"
#include "Explicit_TI_RK.hh"

using namespace Amanzi;

// ODE: y' = y
class fn1 : public Explicit_TI::fnBase<Vector_type> {
 public:
  void FunctionalTimeDerivative(const double t, const Vector_type& y, Vector_type& y_new)
  {
    y_new.assign(y);
  }
};


class fn2 : public Explicit_TI::fnBase<Vector_type> {
 public:
  void FunctionalTimeDerivative(const double t, const Vector_type& y, Vector_type& y_new)
  {
    y_new.putScalar(t * t);
  }
};


TEST(Explicit_RK_Euler)
{
  std::cout << "Test: Explicit_RK_Euler" << std::endl;
  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::forward_euler;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t = 0.0;
  // time step
  double h = 0.1;

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), 2.0 * h);
  }
}


TEST(Explicit_RK_Heun)
{
  std::cout << "Test: Explicit_RK_Heun" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::heun_euler;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time and time step
  double t(0.0), h(0.1);

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 2));
  }
}


TEST(Explicit_RK_Midpoint)
{
  std::cout << "Test: Explicit_RK_Midpoint" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::midpoint;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t = 0.0;
  // time step
  double h = 0.1;

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 2));
  }
}


TEST(Explicit_RK_Ralston)
{
  std::cout << "Test: Explicit_RK_Rapson" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::ralston;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t = 0.0;
  // time step
  double h = 0.1;

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 2));
  }
}


TEST(Explicit_TVD_RK3)
{
  std::cout << "Test: Explicit_TVD_RK3" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::tvd_3rd_order;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t(0.0), h(0.1);

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 3));
  }
}


TEST(Explicit_TVD_RK3_Exact)
{
  std::cout << "Test: Explicit_TVD_RK3 (exactness)" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn2 f;
  auto method = Explicit_TI::tvd_3rd_order;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(0.0);

  // initial time
  double t(0.0), h(0.1);

  // integrate
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(0.0, yv(0, 0) - t * t * t / 3, 1e-15);
  }
}


TEST(Explicit_RK_Kutta3D)
{
  std::cout << "Test: Explicit_RK_Kutta3D" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::kutta_3rd_order;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t = 0.0;
  // time step
  double h = 0.1;

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 3));
  }
}


TEST(Explicit_RK_UserDefined)
{
  std::cout << "Test: Explicit_RK_UserDefined" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  int order = 2;
  Teuchos::SerialDenseMatrix<int, double> a(order, order);
  std::vector<double> b(order);
  std::vector<double> c(order);

  a(1, 0) = 1.0;

  b[0] = 0.5;
  b[1] = 0.5;

  c[0] = 0.0;
  c[1] = 1.0;

  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, order, a, b, c, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t = 0.0;
  // time step
  double h = 0.1;

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 2));
  }
}


TEST(Explicit_RK_RK4)
{
  std::cout << "Test: Explicit_RK_RK4" << std::endl;

  auto comm = getDefaultComm();
  auto map = Teuchos::rcp(new Map_type(1, 1, 0, comm));
  Vector_type y(map), y_new(map);

  fn1 f;
  auto method = Explicit_TI::runge_kutta_4th_order;
  Explicit_TI::RK<Vector_type> explicit_time_integrator(f, method, y);

  // initial value
  y.putScalar(1.0);

  // initial time
  double t = 0.0;
  // time step
  double h = 0.1;

  // integrate to t=1.0
  do {
    explicit_time_integrator.TimeStep(t, h, y, y_new);
    t = t + h;
    y = y_new;
  } while (t < 1.0);

  {
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadOnly);
    CHECK_CLOSE(yv(0, 0), exp(t), pow(h, 4));
  }
}
