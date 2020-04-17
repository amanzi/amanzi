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

#include "AmanziComm.hh"
#include "FunctionConstant.hh"
#include "FunctionFactory.hh"
#include "FunctionLinear.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FunctionPolynomial.hh"
#include "FunctionSeparable.hh"
#include "errors.hh"
#include "CompositeVectorSpace.hh"
#include "MeshFunction.hh"

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
  Kokkos::initialize(argc, argv);
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}


TEST_FIXTURE(reference_mesh, basic_patch)
{
  Teuchos::RCP<MultiFunction> f1 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));

  PatchSpace ps(mesh, true, "LEFT", AmanziMesh::FACE, 1, 1);
  Patch p(ps);

  computeMeshFunction(*f1, 0.0, p);
  CHECK_EQUAL(4, p.size());
  CHECK_EQUAL(1.0, p.data(0,0));
  CHECK_EQUAL(1.0, p.data(3,0));
}


TEST_FIXTURE(reference_mesh, values1)
{
  Teuchos::RCP<MultiFunction> f1 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));
  Teuchos::RCP<MultiFunction> f2 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(2.0))));
  Teuchos::RCP<MultiFunction> f3 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(3.0))));

  MultiPatchSpace mps(mesh, 1, false);
  mps.AddPatch("RIGHT", AmanziMesh::FACE, 1);
  mps.AddPatch("FRONT", AmanziMesh::FACE, 1);
  mps.AddPatch("BACK", AmanziMesh::FACE, 1);
  MultiPatch mp(mps);

  auto funcs = std::vector<Teuchos::RCP<const MultiFunction>>{f1,f2,f3};

  computeMeshFunction(funcs, 0.0, mp);

  CHECK_EQUAL(1.0, mp[0].data(3,0));
  CHECK_EQUAL(2.0, mp[1].data(3,0));
  CHECK_EQUAL(3.0, mp[2].data(3,0));
}


TEST_FIXTURE(reference_mesh, into_vector)
{
  Teuchos::RCP<MultiFunction> f1 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));
  Teuchos::RCP<MultiFunction> f2 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(2.0))));
  Teuchos::RCP<MultiFunction> f3 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(3.0))));

  MultiPatchSpace mps(mesh, 1, false);
  mps.AddPatch("RIGHT", AmanziMesh::FACE, 1);
  mps.AddPatch("FRONT", AmanziMesh::FACE, 1);
  mps.AddPatch("BACK", AmanziMesh::FACE, 1);

  auto funcs = std::vector<Teuchos::RCP<const MultiFunction>>{f1,f2,f3};
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)
      ->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  CompositeVector cv(cvs.CreateSpace());
  cv.putScalar(-1.0);  
  computeMeshFunction(funcs, 0.0, mps, cv);

  {  
    auto cv_v = cv.ViewComponent<MirrorHost>("face", true);
    Entity_ID_List face_list;

    mesh->get_set_entities("RIGHT", FACE, Parallel_type::ALL, face_list);
    for (int i = 0; i < face_list.size(); ++i)
      CHECK_EQUAL(1.0, cv_v(face_list[i], 0));

    mesh->get_set_entities("FRONT", FACE, Parallel_type::ALL, face_list);
    for (int i = 0; i < face_list.size(); ++i)
      CHECK_EQUAL(2.0, cv_v(face_list[i], 0));

    mesh->get_set_entities("BACK", FACE, Parallel_type::ALL, face_list);
    for (int i = 0; i < face_list.size(); ++i)
      CHECK_EQUAL(3.0, cv_v(face_list[i], 0));
  }  
}
