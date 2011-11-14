#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "flow-bc-factory.hh"
#include "boundary-function.hh"
#include "function-factory.hh"
#include "errors.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests ();
}

struct bits_and_pieces
{
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;
  
  enum Side { LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP };
  
  bits_and_pieces()
  {
    comm = new Epetra_MpiComm(MPI_COMM_WORLD);
    MeshFactory mesh_fact(*comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);
  }
};

TEST_FIXTURE(bits_and_pieces, pressure_empty)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory mesh_fact(comm);
  Teuchos::RCP<Mesh> mesh(mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction *bc = bc_fact.CreatePressure();
  bc->Compute(0.0);
  CHECK(bc->end() == bc->begin());
}

TEST_FIXTURE(bits_and_pieces, pressure)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  Teuchos::ParameterList &dir = params->sublist("pressure");
  Teuchos::Array<int> foo_reg(2); foo_reg[0] = LEFT; foo_reg[1] = RIGHT;
  Teuchos::Array<int> bar_reg(1,TOP);
  dir.sublist("foo").set("regions",foo_reg).sublist("boundary pressure").sublist("function-constant").set("value",1.0);
  dir.sublist("bar").set("regions",bar_reg).sublist("boundary pressure").sublist("function-constant").set("value",2.0);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction *bc = bc_fact.CreatePressure();
  bc->Compute(0.0);
  CHECK_EQUAL(12, bc->size());
}

SUITE(pressure_bad_param) {
  TEST_FIXTURE(bits_and_pieces, pressure_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->set("pressure",0); // wrong -- this should be a sublist
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(),Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, spec_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->sublist("pressure").set("fubar", 0); // wrong -- expecting only sublists
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_region)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    foo.sublist("boundary pressure").sublist("function-constant").set("value",0.0);
    // wrong - missing Regions parameter
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
    foo.set("regions",0.0); // wrong -- type should be Array<int>
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_function)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    Teuchos::Array<int> foo_reg(2); foo_reg[0] = LEFT; foo_reg[1] = RIGHT;
    foo.set("regions",foo_reg);
    // wrong - missing boundary pressure list
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
    foo.set("boundary pressure",0); // wrong - not a sublist
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
    foo.remove("boundary pressure");
    foo.sublist("boundary pressure").sublist("function-constant"); // incomplete
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
  }
}

TEST_FIXTURE(bits_and_pieces, mass_flux_empty)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory mesh_fact(comm);
  Teuchos::RCP<Mesh> mesh(mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction *bc = bc_fact.CreatePressure();
  bc->Compute(0.0);
  CHECK(bc->end() == bc->begin());
}

TEST_FIXTURE(bits_and_pieces, mass_flux)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  Teuchos::ParameterList &dir = params->sublist("pressure");
  Teuchos::Array<int> foo_reg(2); foo_reg[0] = LEFT; foo_reg[1] = RIGHT;
  Teuchos::Array<int> bar_reg(1,TOP);
  dir.sublist("foo").set("regions",foo_reg).sublist("boundary pressure").sublist("function-constant").set("value",1.0);
  dir.sublist("bar").set("regions",bar_reg).sublist("boundary pressure").sublist("function-constant").set("value",2.0);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction *bc = bc_fact.CreatePressure();
  bc->Compute(0.0);
  CHECK_EQUAL(12, bc->size());
}

SUITE(mass_flux_bad_param) {
  TEST_FIXTURE(bits_and_pieces, pressure_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->set("pressure",0); // wrong -- this should be a sublist
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(),Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, spec_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->sublist("pressure").set("fubar", 0); // wrong -- expecting only sublists
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_region)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    foo.sublist("boundary pressure").sublist("function-constant").set("value",0.0);
    // wrong - missing Regions parameter
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
    foo.set("regions",0.0); // wrong -- type should be Array<int>
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_function)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("pressure").sublist("foo");
    Teuchos::Array<int> foo_reg(2); foo_reg[0] = LEFT; foo_reg[1] = RIGHT;
    foo.set("regions",foo_reg);
    // wrong - missing boundary pressure list
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
    foo.set("boundary pressure",0); // wrong - not a sublist
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
    foo.remove("boundary pressure");
    foo.sublist("boundary pressure").sublist("function-constant"); // incomplete
    //BoundaryFunction *bc = bc_fact.CreatePressure();
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreatePressure(), Errors::Message);
  }
}

TEST_FIXTURE(bits_and_pieces, static_head_empty)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0, 1.0, 1.0);
  bc->Compute(0.0);
  CHECK(bc->end() == bc->begin());
}

TEST_FIXTURE(bits_and_pieces, static_head)
{
  Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
  Teuchos::ParameterList &dir = params->sublist("static head");
  Teuchos::Array<int> foo_reg(2); foo_reg[0] = LEFT; foo_reg[1] = RIGHT;
  Teuchos::Array<int> bar_reg(1,TOP);
  dir.sublist("foo").set("regions",foo_reg).sublist("water table elevation").sublist("function-constant").set("value",1.0);
  dir.sublist("bar").set("regions",bar_reg).sublist("water table elevation").sublist("function-constant").set("value",2.0);
  FlowBCFactory bc_fact(mesh, params);
  BoundaryFunction *bc0 = bc_fact.CreateStaticHead(0.0, 1.0, 2.0);
  BoundaryFunction *bc1 = bc_fact.CreateStaticHead(1.0, 1.0, 2.0);
  BoundaryFunction *bc2 = bc_fact.CreateStaticHead(0.0, 2.0, 1.0);
  BoundaryFunction *bc3 = bc_fact.CreateStaticHead(0.0, 2.0, 2.0);
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
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0, 1.0, 1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0, 1.0, 1.0),Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, spec_not_list)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    params->sublist("static head").set("fubar", 0); // wrong -- expecting only sublists
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_region)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("static head").sublist("foo");
    foo.sublist("water table elevation").sublist("function-constant").set("value",0.0);
    // wrong - missing Regions parameter
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0), Errors::Message);
    foo.set("regions",0.0); // wrong -- type should be Array<int>
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0), Errors::Message);
  }
  TEST_FIXTURE(bits_and_pieces, bad_function)
  {
    Teuchos::RCP<Teuchos::ParameterList> params(new Teuchos::ParameterList);
    Teuchos::ParameterList &foo = params->sublist("static head").sublist("foo");
    Teuchos::Array<int> foo_reg(2); foo_reg[0] = LEFT; foo_reg[1] = RIGHT;
    foo.set("regions",foo_reg);
    // wrong - missing water table elevation list
    FlowBCFactory bc_fact(mesh, params);
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0), Errors::Message);
    foo.set("water table elevation",0); // wrong - not a sublist
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0), Errors::Message);
    foo.remove("water table elevation");
    foo.sublist("water table elevation").sublist("function-constant"); // incomplete
    //BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0);
    CHECK_THROW(BoundaryFunction *bc = bc_fact.CreateStaticHead(1.0,1.0,1.0), Errors::Message);
  }
}

