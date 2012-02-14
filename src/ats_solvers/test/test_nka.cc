#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"
#include "Mesh_STK.hh"
#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "NonlinarKrylovAccelerator.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
  // data structures for testing
  struct test_data {
    Epetra_MpiComm *comm;
    Teuchos::RCP<AmanziMesh::Mesh> mesh;

    Teuchos::RCP<CompositeVector> x;
    Teuchos::RCP<TreeVector> x_tree;
    Teuchos::RCP<CompositeVector> f;
    Teuchos::RCP<TreeVector> f_tree;
    Teuchos::RCP<CompositeVector> dx;
    Teuchos::RCP<TreeVector> dx_tree;

    Epetra_Vector* fe;
    Epetra_Vector* xe;

    test_data() {
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      AmanziMesh::MeshFactory mesh_fact(*comm);
      mesh = mesh_fact(0.0, 0.0, 0.0, 3.0, 1.0, 1.0, 3, 1, 1);

      // non-ghosted x
      x = Teuchos::rcp(new CompositeVector(mesh, "cell", AmanziMesh::CELL, 1, false));
      x_tree = Teuchos::rcp(new TreeVector("x"));
      x_tree->set_data(x);

      // f and dx
      f_tree = Teuchos::rcp(new TreeVector("f", *x_tree));
      f = f_tree->data();
      dx_tree = Teuchos::rcp(new TreeVector("dx", *x_tree));
      dx = dx_tree->data();

      // just the vector for easier computing
      fe = (*f->ViewComponent("cell"))(0);
      xe = (*x->ViewComponent("cell"))(0);
    }
    ~test_data() {
      delete comm;
    }
  };

  TEST_FIXTURE(test_data, NKA_NONLINEAR) {

    std::cout << "NKA_NONLINEAR..." << std::endl;

    Amanzi::NonlinearKrylovAccelerator FPA(10,0.0,*x_tree);
    FPA.nka_restart();
    FPA.nka_relax();

    // initial value
    x->PutScalar(0.0);

    double norm(0.0);
    int nka_iterations(0);
    do {
      nka_iterations++;

      // function evaluation f <-- f(x^n)
      (*fe)[0] = (*xe)[0] - (cos((*xe)[0]) - sin((*xe)[1])) / 3.0;
      (*fe)[1] = (*xe)[1] - (cos((*xe)[0]) - 2.0*sin((*xe)[1])) / 3.0;
      (*fe)[2] = (*xe)[2] - (cos((*xe)[0]) - 3.0*sin((*xe)[2])) / 3.0;

      // compute the NKA correction  dx <-- NKA(f)
      FPA.nka_correction(*f_tree, *dx_tree);

      // do the fixed point step   x^{n+1} <-- x^n - dx
      x_tree->Update(-1.0, *dx_tree, 1.0);

      dx_tree->NormInf(&norm);

      std::cout << std::setprecision(15) << std::scientific;
      std::cout << "iterate = " << nka_iterations;
      std::cout << ", error = " << norm << std::endl;
    } while (norm > 1e-14);

    std::cout << "final error = " << norm << std::endl;
    CHECK_EQUAL(nka_iterations,9);
    std::cout << "NKA_NONLINEAR... DONE." << std::endl;
  }

  TEST_FIXTURE(test_data, NKA_LINEAR) {
    // This is a linear problem, NKA solves this in N steps where
    // N is the number of unknowns.

    std::cout << "NKA_LINEAR..." << std::endl;

    Amanzi::NonlinearKrylovAccelerator FPA(3,0.01,*x_tree);
    FPA.nka_restart();

    (*xe)[0] = 1.0;
    (*xe)[1] = 2.0;
    (*xe)[2] = -1.0;

    double norm(0.0);
    int nka_iterations(0);
    do {
      nka_iterations++;
      // function evaluation f <-- f(x^n)

      (*fe)[0] = 1.0 - (3.0*(*xe)[0]      -(*xe)[1]      -(*xe)[2]);
      (*fe)[1] = 2.0 - (   -(*xe)[0] + 3.0*(*xe)[1]      -(*xe)[2]);
      (*fe)[2] = -1.0 -(   -(*xe)[0]      -(*xe)[1] + 3.0*(*xe)[2]);

      // compute the NKA correction  dx <-- NKA(f)
      FPA.nka_correction(*f_tree, *dx_tree);

      // do the fixed point step   x^{n+1} <-- x^n - dx
      x_tree->Update(-1.0, *dx_tree, 1.0);

      dx_tree->NormInf(&norm);

      std::cout << std::setprecision(15) << std::scientific;

      std::cout << "iterate = " << nka_iterations;
      std::cout << ", error = " << norm << std::endl;
    } while (norm > 1e-14);

    std::cout << "final error = " << norm << std::endl;
    CHECK_EQUAL(nka_iterations,3);
    std::cout << "NKA_LINEAR... DONE." << std::endl;
  }
}
