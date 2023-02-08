/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "BoundaryFunction.hh"
#include "FunctionConstant.hh"
#include "FunctionFactory.hh"
#include "FunctionLinear.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "FunctionPolynomial.hh"
#include "FunctionSeparable.hh"
#include "errors.hh"

#include "VerboseObject_objs.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto result = UnitTest::RunAllTests();
  Kokkos::finalize();
  return result;
}


struct reference_mesh {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  std::string LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP, INVALID;

  reference_mesh()
  {
    LEFT = "LEFT";
    RIGHT = "RIGHT";
    FRONT = "FRONT";
    BACK = "BACK";
    BOTTOM = "BOTTOM";
    TOP = "TOP";
    INVALID = "INVALID";

    comm = getDefaultComm();

    // Brick domain corners and outward normals to sides
    Teuchos::Array<double> corner_min(Teuchos::tuple(0.0, 0.0, 0.0));
    Teuchos::Array<double> corner_max(Teuchos::tuple(4.0, 4.0, 4.0));
    Teuchos::Array<double> left(Teuchos::tuple(-1.0, 0.0, 0.0));
    Teuchos::Array<double> right(Teuchos::tuple(1.0, 0.0, 0.0));
    Teuchos::Array<double> front(Teuchos::tuple(0.0, -1.0, 0.0));
    Teuchos::Array<double> back(Teuchos::tuple(0.0, 1.0, 0.0));
    Teuchos::Array<double> bottom(Teuchos::tuple(0.0, 0.0, -1.0));
    Teuchos::Array<double> top(Teuchos::tuple(0.0, 0.0, 1.0));

    // Create the geometric model
    Teuchos::ParameterList regions;
    regions.sublist("LEFT").sublist("region: plane").set("point", corner_min).set("normal", left);
    regions.sublist("FRONT").sublist("region: plane").set("point", corner_min).set("normal", front);
    regions.sublist("BOTTOM")
      .sublist("region: plane")
      .set("point", corner_min)
      .set("normal", bottom);
    regions.sublist("RIGHT").sublist("region: plane").set("point", corner_max).set("normal", right);
    regions.sublist("BACK").sublist("region: plane").set("point", corner_max).set("normal", back);
    regions.sublist("TOP").sublist("region: plane").set("point", corner_max).set("normal", top);
    Teuchos::RCP<AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new AmanziGeometry::GeometricModel(3, regions, *comm));
    // Create the mesh
    MeshFactory mesh_fact(comm, gm);
    mesh = mesh_fact.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
  }
};


TEST_FIXTURE(reference_mesh, empty)
{
  BoundaryFunction bf(mesh);
  CHECK_EQUAL(0, bf.size());
}


TEST_FIXTURE(reference_mesh, basic)
{
  BoundaryFunction bf(mesh);
  Teuchos::RCP<MultiFunction> f1 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));
  // Add a definition for a single side.
  bf.Define(RIGHT, f1);
  // Add a definition for a couple other sides
  std::vector<std::string> reg(2);
  reg[0] = FRONT;
  reg[1] = BACK;
  bf.Define(reg, f1);
  // Add a definition for yet more sides, but with duplicates.
  reg[0] = BOTTOM;
  reg[1] = BOTTOM;
  bf.Define(reg, f1);
  bf.Finalize();
  CHECK_EQUAL(16, bf.size());
}


TEST_FIXTURE(reference_mesh, values1)
{
  BoundaryFunction bf(mesh);
  Teuchos::RCP<MultiFunction> f1 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));
  Teuchos::RCP<MultiFunction> f2 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(2.0))));
  Teuchos::RCP<MultiFunction> f3 =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(3.0))));

  bf.Define(RIGHT, std::move(f1));
  bf.Define(FRONT, std::move(f2));
  bf.Define(BACK, std::move(f3));
  bf.Finalize();

  CHECK_EQUAL(12, bf.size());
  bf.Compute(0.0);

  auto face_list = mesh->getSetEntities(RIGHT, FACE, Parallel_kind::ALL);
  for (auto f = face_list.begin(); f != face_list.end(); ++f) CHECK_EQUAL(1.0, bf.find(*f)->second);

  face_list = mesh->getSetEntities(FRONT, FACE, Parallel_kind::ALL);
  for (auto f = face_list.begin(); f != face_list.end(); ++f) CHECK_EQUAL(2.0, bf.find(*f)->second);

  face_list = mesh->getSetEntities(BACK, FACE, Parallel_kind::ALL);
  for (auto f = face_list.begin(); f != face_list.end(); ++f) CHECK_EQUAL(3.0, bf.find(*f)->second);
}


TEST_FIXTURE(reference_mesh, values2)
{
  // Create the function f(t,x,y,z) = t * (x + 2y + 3z)
  AmanziMesh::Double_List c(1, 1.0);
  std::vector<int> ip(1, 1);
  std::unique_ptr<Function> f1(new FunctionPolynomial(c, ip));
  double g[3] = { 1.0, 2.0, 3.0 };
  AmanziMesh::Double_List grad(g, g + 3);
  std::unique_ptr<Function> f2(new FunctionLinear(0.0, grad));

  // Create the boundary function
  Teuchos::RCP<MultiFunction> f3 = Teuchos::rcp(
    new MultiFunction(Teuchos::rcp(new FunctionSeparable(std::move(f1), std::move(f2)))));
  BoundaryFunction bf(mesh);
  std::vector<std::string> regions(2);
  regions[0] = RIGHT;
  regions[1] = BACK;
  bf.Define(regions, f3);

  // Check values at t=1
  bf.Compute(1.0);
  for (BoundaryFunction::Iterator i = bf.begin(); i != bf.end(); ++i) {
    AmanziGeometry::Point p = mesh->getFaceCentroid(i->first);
    CHECK_EQUAL(p.x() + 2 * p.y() + 3 * p.z(), i->second);
  }

  // Check values at t=2
  bf.Compute(2.0);
  for (BoundaryFunction::Iterator i = bf.begin(); i != bf.end(); ++i) {
    AmanziGeometry::Point p = mesh->getFaceCentroid(i->first);
    CHECK_EQUAL(2 * (p.x() + 2 * p.y() + 3 * p.z()), i->second);
  }
}


TEST_FIXTURE(reference_mesh, bad_input)
{
  Teuchos::RCP<MultiFunction> f =
    Teuchos::rcp(new MultiFunction(Teuchos::rcp(new FunctionConstant(1.0))));
  BoundaryFunction bf(mesh);
  //bf.Define(INVALID, f); // no such face set
  CHECK_THROW(bf.Define(INVALID, f), Errors::Message);

  bf.Define(RIGHT, f);
  //bf.Define(RIGHT, f); // overlapping definition
  CHECK_THROW(bf.Define(RIGHT, f), Errors::Message);
}
