/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   Amanzi

   Unit tests for the composite vector.
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
    comm = Amanzi::getDefaultComm();
    MeshFactory mesh_factory(comm);
    mesh_factory.set_preference(Preference({ Framework::MSTK }));
    mesh = mesh_factory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    CompositeVectorSpace cvs;
    cvs.SetMesh(mesh)
      ->SetGhosted(true)
      ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1)
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 2);
    x = Teuchos::rcp(new CompositeVector(cvs));
  }
  ~test_cv() {}
};

SUITE(COMPOSITE_VECTOR)
{
  TEST_FIXTURE(test_cv, CVAccessTiming)
  {
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

    Teuchos::RCP<Epetra_MultiVector> mv = x->viewComponent("cell", false);
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
