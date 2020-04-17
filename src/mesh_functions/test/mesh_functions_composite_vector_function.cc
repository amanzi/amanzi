/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "MultiFunction.hh"
#include "FunctionConstant.hh"
#include "FunctionLinear.hh"
#include "CompositeVectorFunction.hh"
#include "errors.hh"
#include "CompositeVectorSpace.hh"

#include "VerboseObject_objs.hh"

#include "test/reference_mesh.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}


TEST_FIXTURE(reference_mesh, cv_function)
{
  // make the mesh function
  Kokkos::View<double*> grad("grad", 4);
  grad(0) = 0.0;
  grad(1) = 0.0;
  grad(2) = 0.0;
  grad(3) = 0.0;

  Teuchos::RCP<const Function> linear_func =
    Teuchos::rcp(new FunctionLinear(1.0, grad));
  std::vector<Teuchos::RCP<const Function>> linear_funcs(1, linear_func);
  Teuchos::RCP<MultiFunction> f1 = Teuchos::rcp(new MultiFunction(linear_funcs));

  MultiPatchSpace mps(mesh, true, 1);
  mps.AddPatch("DOMAIN", AmanziMesh::CELL, 1);
  mps.AddPatch("DOMAIN", AmanziMesh::FACE, 1);
  mps.AddPatch("DOMAIN", AmanziMesh::NODE, 1);

  auto funcs = std::vector<Teuchos::RCP<const MultiFunction>>{ f1, f1, f1 };

  CompositeVectorFunction cv_func(mps, funcs);
  
  // make the CV
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("face", AmanziMesh::FACE, 1);

  auto cv = cvs.Create();
  cv->putScalar(0.0);

  // apply the function to the vector
  cv_func.Compute(0.0, *cv);

  // Check
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  {
    auto cv_c = cv->ViewComponent<MirrorHost>("cell", false);
    for (int c=0; c!=ncells; ++c) {
      CHECK_CLOSE(1.0, cv_c(c, 0), 0.0000001);
    }
  }

  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  {
    auto cv_c = cv->ViewComponent<MirrorHost>("face", false);
    for (int c=0; c!=nfaces; ++c) {
      CHECK_CLOSE(1.0, cv_c(c, 0), 0.0000001);
    }
  }
}
