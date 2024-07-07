/*
  Copyright 2010-202x held jointly by participating institutions.
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
#include "Teuchos_ParameterXMLFileReader.hpp"

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

#include "test/reference_mesh.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;


TEST_FIXTURE(reference_mesh, COMPOSITE_VECTOR_FUNCTION)
{
  // make the mesh function
  Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 4);
  grad(0) = 0.0;
  grad(1) = 0.0;
  grad(2) = 0.0;
  grad(3) = 0.0;

  Teuchos::RCP<const Function> linear_func = Teuchos::rcp(new FunctionLinear(1.0, grad));
  std::vector<Teuchos::RCP<const Function>> linear_funcs(1, linear_func);
  Teuchos::RCP<MultiFunction> f1 = Teuchos::rcp(new MultiFunction(linear_funcs));

  CompositeVectorFunction cv_func(mesh);
  cv_func.addSpec("cell", AmanziMesh::Entity_kind::CELL, 1, "DOMAIN", f1);
  cv_func.addSpec("face", AmanziMesh::Entity_kind::FACE, 1, "DOMAIN", f1);

  // make the CV
  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh)
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("face", AmanziMesh::Entity_kind::FACE, 1);

  auto cv = cvs.Create();
  cv->putScalar(0.0);

  // apply the function to the vector
  cv_func.Compute(0.0, *cv);

  // Check
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  {
    auto cv_c = cv->viewComponent<MemSpace_kind::HOST>("cell", false);
    for (int c = 0; c != ncells; ++c) { CHECK_CLOSE(1.0, cv_c(c, 0), 0.0000001); }
  }

  int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  {
    auto cv_c = cv->viewComponent<MemSpace_kind::HOST>("face", false);
    for (int c = 0; c != nfaces; ++c) { CHECK_CLOSE(1.0, cv_c(c, 0), 0.0000001); }
  }
}


TEST_FIXTURE(reference_mesh, COMPOSITE_VECTOR_FUNCTION_PLIST)
{
  std::string xmlFileName = "test/mesh_functions_composite_vector_function.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = xmlreader.getParameters();

  CompositeVectorFunction cv_func(plist, mesh);
  CompositeVector cv(cv_func.createCVS(false)->CreateSpace());
  cv_func.Compute(0.0, cv);

  // Check
  int ncells = mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  {
    auto cv_c = cv.viewComponent<MemSpace_kind::HOST>("cell", false);
    CHECK_CLOSE(0.5, cv_c(0, 0), 0.0000001);
    CHECK_CLOSE(2.0, cv_c(ncells - 1, 0), 0.0000001);
  }

  int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  {
    auto cv_c = cv.viewComponent<MemSpace_kind::HOST>("face", false);
    for (int c = 0; c != nfaces; ++c) { CHECK_CLOSE(3.0, cv_c(c, 0), 0.0000001); }
  }
}
