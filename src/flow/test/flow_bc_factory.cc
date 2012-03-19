#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "Flow_BC_Factory.hpp"
#include "boundary-function.hh"
#include "function-factory.hh"
#include "errors.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::AmanziFlow;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests ();
}


struct bits_and_pieces
{
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;
  GeometricModel *gm;
  
  enum Side {LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP};
  
  bits_and_pieces()
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
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
    regions.sublist("LEFT").sublist("Region: Plane").
        set("Location",corner_min).set("Direction",left);
    regions.sublist("FRONT").sublist("Region: Plane").
        set("Location",corner_min).set("Direction",front);
    regions.sublist("BOTTOM").sublist("Region: Plane").
        set("Location",corner_min).set("Direction",bottom);
    regions.sublist("RIGHT").sublist("Region: Plane").
        set("Location",corner_max).set("Direction",right);
    regions.sublist("BACK").sublist("Region: Plane").
        set("Location",corner_max).set("Direction",back);
    regions.sublist("TOP").sublist("Region: Plane").
        set("Location",corner_max).set("Direction",top);
    gm = new GeometricModel(3, regions, comm);
    // Create the mesh
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2, gm);
  }
};


TEST_FIXTURE(bits_and_pieces, pressure_empty)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory mesh_fact(&comm);
  Teuchos::RCP<Mesh> mesh(mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction* bc = bc_fact.createPressure();
  bc->Compute(0.0);
  CHECK(bc->end() == bc->begin());
  delete bc;
}


TEST_FIXTURE(bits_and_pieces, pressure)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  Teuchos::ParameterList &dir = params->sublist("pressure");
  Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"),std::string("RIGHT")));
  Teuchos::Array<std::string> bar_reg(Teuchos::tuple(std::string("TOP")));
  dir.sublist("foo").set("regions",foo_reg).sublist("boundary pressure").sublist("function-constant").set("value",1.0);
  dir.sublist("bar").set("regions",bar_reg).sublist("boundary pressure").sublist("function-constant").set("value",2.0);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction* bc = bc_fact.createPressure();
  bc->Compute(0.0);
  CHECK_EQUAL(12, bc->size());
  delete bc;
}


SUITE(pressure_bad_param) {
  TEST_FIXTURE(bits_and_pieces, pressure_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->set("pressure", 0); // wrong -- this should be a sublist
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }

  TEST_FIXTURE(bits_and_pieces, spec_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->sublist("pressure").set("fubar", 0); // wrong -- expecting only sublists
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }

  TEST_FIXTURE(bits_and_pieces, bad_region)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    foo.sublist("boundary pressure").sublist("function-constant").set("value",0.0);
    // wrong - missing Regions parameter
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
    foo.set("regions",0.0); // wrong -- type should be Array<string>
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }

  TEST_FIXTURE(bits_and_pieces, bad_function)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"),std::string("RIGHT")));
    foo.set("regions",foo_reg);
    // wrong - missing boundary pressure list
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
    foo.set("boundary pressure",0); // wrong - not a sublist
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
    foo.remove("boundary pressure");
    foo.sublist("boundary pressure").sublist("function-constant"); // incomplete
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }
}

TEST_FIXTURE(bits_and_pieces, mass_flux_empty)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory mesh_fact(&comm);
  Teuchos::RCP<Mesh> mesh(mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction* bc = bc_fact.createPressure();
  bc->Compute(0.0);
  CHECK(bc->end() == bc->begin());
}

TEST_FIXTURE(bits_and_pieces, mass_flux)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  Teuchos::ParameterList &dir = params->sublist("pressure");
  Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"),std::string("RIGHT")));
  Teuchos::Array<std::string> bar_reg(Teuchos::tuple(std::string("TOP")));
  dir.sublist("foo").set("regions",foo_reg).sublist("boundary pressure").sublist("function-constant").set("value",1.0);
  dir.sublist("bar").set("regions",bar_reg).sublist("boundary pressure").sublist("function-constant").set("value",2.0);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction* bc = bc_fact.createPressure();
  bc->Compute(0.0);
  CHECK_EQUAL(12, bc->size());
}

SUITE(mass_flux_bad_param) {
  TEST_FIXTURE(bits_and_pieces, pressure_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->set("pressure",0); // wrong -- this should be a sublist
    FlowBCFactory bc_fact(mesh, params);
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(),Errors::Message);
  }

  TEST_FIXTURE(bits_and_pieces, spec_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->sublist("pressure").set("fubar", 0); // wrong -- expecting only sublists
    FlowBCFactory bc_fact(mesh, params);
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }

  TEST_FIXTURE(bits_and_pieces, bad_region)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    foo.sublist("boundary pressure").sublist("function-constant").set("value",0.0);
    // wrong - missing Regions parameter
    FlowBCFactory bc_fact(mesh, params);
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
    foo.set("regions",0.0); // wrong -- type should be Array<string>
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }

  TEST_FIXTURE(bits_and_pieces, bad_function)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"),std::string("RIGHT")));
    foo.set("regions",foo_reg);
    // wrong - missing boundary pressure list
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
    foo.set("boundary pressure",0); // wrong - not a sublist
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
    foo.remove("boundary pressure");
    foo.sublist("boundary pressure").sublist("function-constant"); // incomplete
    //BoundaryFunction* bc = bc_fact.createPressure();
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createPressure(), Errors::Message);
  }
}

TEST_FIXTURE(bits_and_pieces, static_head_empty)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  FlowBCFactory bc_fact(mesh, params);
  AmanziGeometry::Point gravity(0.0, 0.0, -1.0);

  BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, gravity);
  bc->Compute(0.0);
  CHECK(bc->end() == bc->begin());
}

TEST_FIXTURE(bits_and_pieces, static_head)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  Teuchos::ParameterList &dir = params->sublist("static head");
  Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"),std::string("RIGHT")));
  Teuchos::Array<std::string> bar_reg(Teuchos::tuple(std::string("TOP")));
  dir.sublist("foo").set("regions",foo_reg).sublist("water table elevation").sublist("function-constant").set("value",1.0);
  dir.sublist("bar").set("regions",bar_reg).sublist("water table elevation").sublist("function-constant").set("value",2.0);
  FlowBCFactory bc_fact(mesh, params);
  AmanziGeometry::Point g2(0.0, 0.0, -2.0), g1(0.0, 0.0, -1.0);

  BoundaryFunction* bc0 = bc_fact.createStaticHead(0.0, 1.0, g2);
  BoundaryFunction* bc1 = bc_fact.createStaticHead(1.0, 1.0, g2);
  BoundaryFunction* bc2 = bc_fact.createStaticHead(0.0, 2.0, g1);
  BoundaryFunction *bc3 = bc_fact.createStaticHead(0.0, 2.0, g2);
  bc0->Compute(0.0);
  CHECK_EQUAL(12, bc0->size());

  BoundaryFunction::Iterator i, j;
  bc1->Compute(0.0);
  for (i = bc0->begin(), j = bc1->begin(); i != bc0->end(); ++i, ++j) {
    CHECK_EQUAL(1+ i->second, j->second);
  }
  bc2->Compute(0.0);
  for (i = bc0->begin(), j = bc2->begin(); i != bc0->end(); ++i, ++j) {
    CHECK_EQUAL(i->second, j->second);
  }
  bc3->Compute(0.0);
  for (i = bc0->begin(), j = bc3->begin(); i != bc0->end(); ++i, ++j) {
    CHECK_EQUAL(2*i->second, j->second);
  }
}

SUITE(static_head_bad_param) {
  TEST_FIXTURE(bits_and_pieces, static_head_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->set("static head",0); // wrong -- this should be a sublist
    FlowBCFactory bc_fact(mesh, params);
    AmanziGeometry::Point g(0.0, 0.0, -1.0);
    CHECK_THROW(BoundaryFunction*bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, spec_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->sublist("static head").set("fubar", 0); // wrong -- expecting only sublists
    FlowBCFactory bc_fact(mesh, params);
    AmanziGeometry::Point g(0.0, 0.0, -1.0);
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_region)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("static head").sublist("foo");
    foo.sublist("water table elevation").sublist("function-constant").set("value",0.0);
    // wrong - missing Regions parameter
    FlowBCFactory bc_fact(mesh, params);
    AmanziGeometry::Point g(0.0, 0.0, -1.0);
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
    foo.set("regions", 0.0); // wrong -- type should be Array<string>
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_function)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("static head").sublist("foo");
    Teuchos::Array<std::string> foo_reg(Teuchos::tuple(std::string("LEFT"),std::string("RIGHT")));
    foo.set("regions",foo_reg);
    // wrong - missing water table elevation list
    FlowBCFactory bc_fact(mesh, params);
    AmanziGeometry::Point g(0.0, 0.0, -1.0);
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
    foo.set("water table elevation", 0); // wrong - not a sublist
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
    foo.remove("water table elevation");
    foo.sublist("water table elevation").sublist("function-constant"); // incomplete
    CHECK_THROW(BoundaryFunction* bc = bc_fact.createStaticHead(1.0, 1.0, g), Errors::Message);
  }
}

