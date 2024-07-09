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


TEST_FIXTURE(reference_mesh, MESH_FUNCTION)
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

  MeshFunction mf(mesh, AmanziMesh::Entity_kind::CELL);
  mf.addSpec("cell", AmanziMesh::Entity_kind::CELL, 1, "HALF1", f1);
  mf.addSpec("cell", AmanziMesh::Entity_kind::CELL, 1, "HALF2", f1);

  auto mps = mf.createMPS(false);
  MultiPatch<double> mp(mps);
  CHECK_EQUAL(2, mp.size());

  mf.Compute(0.0, mp);
  auto mp_host = Kokkos::create_mirror_view_and_copy(DefaultHostMemorySpace(), mp[0].data);
  CHECK_CLOSE(1.0, mp_host(2,0), 1.e-8);
}


TEST_FIXTURE(reference_mesh, BOUNDARY_CONDITION_PLIST)
{
  std::string xmlFileName = "test/mesh_functions_base.xml";
  Teuchos::ParameterXMLFileReader xmlreader(xmlFileName);
  auto plist = xmlreader.getParameters();

  int OPERATOR_BC_NEUMANN = 2; // Operators::OPERATOR_BC_NEUMANN
  MeshFunction mf(
    plist, mesh, "outward normal flux", AmanziMesh::Entity_kind::FACE, OPERATOR_BC_NEUMANN);

  MultiPatch<double> mp(mf.createMPS(false));
  mf.Compute(0.0, mp);

  auto mp0_host = Kokkos::create_mirror_view_and_copy(DefaultHostMemorySpace(), mp[0].data);
  auto mp1_host = Kokkos::create_mirror_view_and_copy(DefaultHostMemorySpace(), mp[1].data);
  CHECK_CLOSE(0.5, mp0_host(0,0), 1.e-8);
  CHECK_CLOSE(2.0, mp1_host(0,0), 1.e-8);
}
