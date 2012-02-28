#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "Mesh_STK.hh"
#include "MeshFactory.hh"
#include "CompositeVector.hh"
#include "TreeVector.hh"

#include "ExplicitTIfnBase.hh"
#include "ExplicitTIRK.hh"


using namespace Amanzi;

SUITE(TimeIntegrationTests) {
  // data structures
  struct test_data {
    Epetra_MpiComm *comm;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;

    Teuchos::RCP<CompositeVector> y;
    Teuchos::RCP<TreeVector> y_tree;
    Teuchos::RCP<CompositeVector> y_new;
    Teuchos::RCP<TreeVector> y_new_tree;

    Epetra_Vector* ye;
    Epetra_Vector* yne;

    test_data() {
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      AmanziMesh::MeshFactory mesh_fact(comm);
      mesh = mesh_fact(0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 2, 1, 1);

      // non-ghosted y
      y = Teuchos::rcp(new CompositeVector(mesh, "cell", AmanziMesh::CELL, 1, false));
      y_tree = Teuchos::rcp(new TreeVector("y"));
      y_tree->set_data(y);

      // y new
      y_new_tree = Teuchos::rcp(new TreeVector("y_new", *y_tree));
      y_new = y_new_tree->data();

      // just the vector for easier computing
      ye = (*y->ViewComponent("cell"))(0);
      yne = (*y_new->ViewComponent("cell"))(0);
    }
    ~test_data() { delete comm; }
  };

  // ODE: y' = y
  class fn1 : public ExplicitTIfnBase {
  public:
    void fun(const double t, const TreeVector& y, TreeVector& y_new)
    {
      y_new = y;
    }
  };

  TEST_FIXTURE(test_data, Explicit_RK_Euler) {
    cout << "Test: Explicit_RK_Euler" << endl;
    fn1 f;
    ExplicitTIRK::method_t method = ExplicitTIRK::forward_euler;
    ExplicitTIRK explicit_time_integrator(f, method, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t=t+h;
      *y_tree = *y_new_tree;
    } while (t<1.0);
    CHECK_CLOSE((*ye)[0], exp(t), 2.0*h);
  }

  TEST_FIXTURE(test_data, Explicit_RK_Heun) {
    cout << "Test: Explicit_RK_Heun" << endl;
    fn1 f;
    ExplicitTIRK::method_t method = ExplicitTIRK::heun_euler;
    ExplicitTIRK explicit_time_integrator(f, method, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t=t+h;
      *y_tree = *y_new_tree;
    } while (t<1.0);
    CHECK_CLOSE((*ye)[0], exp(t), pow(h,2));
  }

  TEST_FIXTURE(test_data, Explicit_RK_Midpoint) {
    cout << "Test: Explicit_RK_Midpoint" << endl;
    fn1 f;
    ExplicitTIRK::method_t method = ExplicitTIRK::midpoint;
    ExplicitTIRK explicit_time_integrator(f, method, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t=t+h;
      *y_tree = *y_new_tree;
    } while (t<1.0);

    CHECK_CLOSE((*ye)[0], exp(t), pow(h,2));
  }

  TEST_FIXTURE(test_data, Explicit_RK_Ralston) {
    cout << "Test: Explicit_RK_Rapson" << endl;
    fn1 f;
    ExplicitTIRK::method_t method = ExplicitTIRK::ralston;
    ExplicitTIRK explicit_time_integrator(f, method, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t=t+h;
      *y_tree = *y_new_tree;
    } while (t<1.0);
    CHECK_CLOSE((*ye)[0], exp(t), pow(h,2));
  }

  TEST_FIXTURE(test_data, Explicit_RK_Kutta3D) {
    cout << "Test: Explicit_RK_Kutta3D" << endl;
    fn1 f;
    ExplicitTIRK::method_t method = ExplicitTIRK::kutta_3rd_order;
    ExplicitTIRK explicit_time_integrator(f, method, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t=t+h;
      *y_tree = *y_new_tree;
    } while (t<1.0);
    CHECK_CLOSE((*ye)[0], exp(t), pow(h,3));
  }

  TEST_FIXTURE(test_data, Explicit_RK_UserDefined) {
    cout << "Test: Explicit_RK_UserDefined" << endl;
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

    ExplicitTIRK explicit_time_integrator(f, order, a, b, c, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t = t + h;
      *y_tree = *y_new_tree;
    } while (t<1.0);

    CHECK_CLOSE((*ye)[0], exp(t), pow(h,2));
  }

  TEST_FIXTURE(test_data, Explicit_RK_RK4) {
    cout << "Test: Explicit_RK_RK4" << endl;
    fn1 f;
    ExplicitTIRK::method_t method = ExplicitTIRK::runge_kutta_4th_order;
    ExplicitTIRK explicit_time_integrator(f, method, *y_tree);

    // initial value
    y_tree->PutScalar(1.0);

    // initial time
    double t=0.0;
    // time step
    double h=.1;

    // integrate to t=1.0
    do
    {
      explicit_time_integrator.step(t, h, *y_tree, *y_new_tree);
      t = t + h;
      *y_tree = *y_new_tree;
    }
    while (t<1.0);
    CHECK_CLOSE((*ye)[0], exp(t), pow(h,4));
  }
}
