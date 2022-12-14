/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   Amanzi

   ------------------------------------------------------------------------- */

#include <vector>

#include "UnitTest++.h"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "MeshFactory.hh"
#include "Mesh_simple.hh"
#include "CompositeVector.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

struct test_cv {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;

  Teuchos::RCP<CompositeVector> x;
  Teuchos::RCP<CompositeVector> x2;

  test_cv()
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    std::vector<Entity_kind> locations(2);
    locations[0] = CELL;
    locations[1] = FACE;

    std::vector<std::string> names(2);
    names[0] = "cell";
    names[1] = "face";

    std::vector<int> num_dofs(2);
    num_dofs[0] = 2;
    num_dofs[1] = 1;

    x = Teuchos::rcp(new CompositeVector(mesh, names, locations, num_dofs, true));
    //    x2 = Teuchos::rcp(new CompositeVector(mesh, CELL, 1, true));
  }
  ~test_cv() {}
};

SUITE(COMPOSITE_VECTOR)
{
  TEST_FIXTURE(test_cv, CVAccessTiming)
  {
    x->CreateData();
    x->PutScalar(2.0);

    Teuchos::RCP<Teuchos::Time> cvtime =
      Teuchos::TimeMonitor::getNewCounter("composite vector access");
    Teuchos::RCP<Teuchos::Time> mvtime = Teuchos::TimeMonitor::getNewCounter("multivector access");

    int ncells = x->size("cell", false);
    double val(0.0);

    if (true) {
      Teuchos::TimeMonitor timer(*cvtime);
      for (int i = 0; i != 10000000; ++i) {
        for (int j = 0; j != ncells; ++j) { val = (*x)("cell", j); }
      }
    }
    std::cout << val;

    Teuchos::RCP<Epetra_MultiVector> mv = x->ViewComponent("cell", false);
    if (true) {
      Teuchos::TimeMonitor timer(*mvtime);
      for (int i = 0; i != 10000000; ++i) {
        val = 0.0;
        for (int j = 0; j != ncells; ++j) { val = (*mv)[0][j]; }
      }
    }
    std::cout << val;

    Teuchos::TimeMonitor::summarize();
  }
}
