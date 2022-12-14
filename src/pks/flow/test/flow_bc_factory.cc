/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson, version 1 (nnc@lanl.gov),
      Konstantin Lipnikov, version 2 (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include "UnitTest++.h"

#include "TestReporterStdout.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "MeshFactory.hh"
#include "MultiFunction.hh"
#include "PK_DomainFunctionFactory.hh"

#include "Flow_PK.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Flow;

struct bits_and_pieces {
  Comm_ptr_type comm;
  Teuchos::RCP<Mesh> mesh;
  Teuchos::RCP<GeometricModel> gm;

  enum Side { LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP };

  bits_and_pieces()
  {
    comm = Amanzi::getDefaultComm();
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
    gm = Teuchos::rcp(new GeometricModel(3, regions, *comm));
    // Create the mesh
    MeshFactory meshfactory(comm, gm);
    mesh = meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
  }
};


TEST_FIXTURE(bits_and_pieces, empty_parameter_list)
{
  Comm_ptr_type comm = Amanzi::getDefaultComm();
  MeshFactory meshfactory(comm, gm);
  // Teuchos::RCP<Mesh> mesh(meshfactory.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));

  PK_DomainFunctionFactory<PK_DomainFunction> bc_fact(mesh, Teuchos::null);

  Teuchos::ParameterList plist;
  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);
}


TEST_FIXTURE(bits_and_pieces, pressure_not_list)
{
  Teuchos::ParameterList plist;
  plist.set("pressure", 0); // wrong -- this should be a sublist

  PK_DomainFunctionFactory<PK_DomainFunction> bc_fact(mesh, Teuchos::null);
  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);
}


TEST_FIXTURE(bits_and_pieces, bad_region)
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& foo = plist.sublist("pressure").sublist("foo");
  foo.sublist("boundary pressure").sublist("function-constant").set("value", 0.0);
  // wrong - missing regions parameter
  PK_DomainFunctionFactory<PK_DomainFunction> bc_fact(mesh, Teuchos::null);

  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);

  foo.set("spatial distribution method", "none");
  foo.set("regions", 0.0); // wrong -- type should be Array<std::string>
  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);
}


TEST_FIXTURE(bits_and_pieces, bad_function)
{
  Teuchos::ParameterList plist;
  Teuchos::ParameterList& foo = plist.sublist("pressure").sublist("foo");
  Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"), std::string("RIGHT")));
  foo.set("spatial distribution method", "none");
  foo.set("regions", foo_reg);
  // wrong - missing boundary pressure list

  PK_DomainFunctionFactory<PK_DomainFunction> bc_fact(mesh, Teuchos::null);
  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);

  foo.set("boundary pressure", 0); // wrong - not a sublist
  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);

  foo.remove("boundary pressure");
  foo.sublist("boundary pressure").sublist("function-constant"); // incomplete
  CHECK_THROW(bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);
}


TEST_FIXTURE(bits_and_pieces, pressure)
{
  Teuchos::ParameterList plist;

  Teuchos::Array<std::string> regs(Teuchos::tuple(std::string("LEFT"), std::string("RIGHT")));
  plist.set("regions", regs)
    .set("spatial distribution method", "none")
    .sublist("boundary pressure")
    .sublist("function-constant")
    .set("value", 1.0);

  PK_DomainFunctionFactory<PK_DomainFunction> bc_fact(mesh, Teuchos::null);
  Teuchos::RCP<PK_DomainFunction> bc =
    bc_fact.Create(plist, "boundary pressure", AmanziMesh::FACE, Teuchos::null);

  bc->Compute(0.0, 0.0);
  CHECK_EQUAL(8, bc->size());
}


TEST_FIXTURE(bits_and_pieces, static_head_empty)
{
  Teuchos::ParameterList plist;
  AmanziGeometry::Point gravity(0.0, 0.0, -1.0);

  PK_DomainFunctionFactory<PK_DomainFunction> bc_fact(mesh, Teuchos::null);
  CHECK_THROW(bc_fact.Create(plist, "static head", AmanziMesh::FACE, Teuchos::null),
              Errors::Message);
}


TEST_FIXTURE(bits_and_pieces, static_head)
{
  Teuchos::ParameterList foo, bar;
  Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"), std::string("RIGHT")));
  Teuchos::Array<std::string> bar_reg(Teuchos::tuple(std::string("TOP")));
  foo.set("regions", foo_reg)
    .set("spatial distribution method", "none")
    .sublist("static head")
    .sublist("function-static-head")
    .set("p0", 0.0)
    .set("density", 1.0)
    .set("gravity", 2.0)
    .set("space dimension", 3)
    .sublist("water table elevation")
    .sublist("function-constant")
    .set("value", 1.0);
  bar.set("regions", bar_reg)
    .set("spatial distribution method", "none")
    .sublist("static head")
    .sublist("function-static-head")
    .set("p0", 1.0)
    .set("density", 1.0)
    .set("gravity", 2.0)
    .set("space dimension", 3)
    .sublist("water table elevation")
    .sublist("function-constant")
    .set("value", 2.0);

  PK_DomainFunctionFactory<FlowBoundaryFunction> bc_fact(mesh, Teuchos::null);

  Teuchos::RCP<FlowBoundaryFunction> bc0 =
    bc_fact.Create(foo, "static head", AmanziMesh::FACE, Teuchos::null);
  Teuchos::RCP<FlowBoundaryFunction> bc1 =
    bc_fact.Create(bar, "static head", AmanziMesh::FACE, Teuchos::null);

  bc0->Compute(0.0, 0.0);
  bc0->ComputeSubmodel(mesh);

  bc1->Compute(0.0, 0.0);
  bc1->ComputeSubmodel(mesh);

  CHECK_EQUAL(8, bc0->size());

  int i;
  double head[8] = { 0.0, -4.0, 0.0, -4.0, 0.0, -4.0, 0.0, -4.0 };
  FlowBoundaryFunction ::Iterator it0, it1;

  for (it0 = bc0->begin(), i = 0; it0 != bc0->end(); ++it0, ++i) {
    CHECK_EQUAL(head[i], it0->second[0]);
  }
  for (it1 = bc1->begin(); it1 != bc1->end(); ++it1) { CHECK_EQUAL(-3.0, it1->second[0]); }
}
