#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "boundary_function.hh"
#include "ConstantFunction.hh"
#include "PolynomialFunction.hh"
#include "LinearFunction.hh"
#include "SeparableFunction.hh"
#include "FunctionFactory.hh"
#include "composite_function.hh"
#include "errors.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::Functions;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests ();
}

struct reference_mesh
{
  Epetra_MpiComm *comm;
  Teuchos::RCP<Mesh> mesh;
  GeometricModel *gm;
  std::string LEFT;
  std::string RIGHT;
  std::string FRONT;
  std::string BACK;
  std::string BOTTOM;
  std::string TOP;
  std::string INVALID;

  reference_mesh()
  {
    LEFT   = "LEFT";
    RIGHT  = "RIGHT";
    FRONT  = "FRONT";
    BACK   = "BACK";
    BOTTOM = "BOTTOM";
    TOP    = "TOP";
    INVALID = "INVALID";

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
    gm = new GeometricModel(3,regions,comm);
    // Create the mesh
    MeshFactory mesh_fact(comm);
    mesh = mesh_fact(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2, gm);
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
  Teuchos::RCP<CompositeFunction> f1 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(1.0))));
  Teuchos::RCP<CompositeFunction> f2 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(2.0))));
  Teuchos::RCP<CompositeFunction> f3 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(3.0))));
  // Add a definition for a single side.
  bf.Define(RIGHT, f1);
  // Add a definition for a couple other sides
  std::vector<std::string> reg(2);
  reg[0] = FRONT; reg[1] = BACK;
  bf.Define(reg, f1);
  // Add a definition for yet more sides, but with duplicates.
  reg[0] = BOTTOM; reg[1] = BOTTOM;
  bf.Define(reg, f1);
  bf.Finalize();
  CHECK_EQUAL(16,bf.size());
}

TEST_FIXTURE(reference_mesh, values1)
{
  BoundaryFunction bf(mesh);
  Teuchos::RCP<CompositeFunction> f1 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(1.0))));
  Teuchos::RCP<CompositeFunction> f2 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(2.0))));
  Teuchos::RCP<CompositeFunction> f3 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(3.0))));
  bf.Define(RIGHT, f1);
  bf.Define(FRONT, f2);
  bf.Define(BACK,  f3);
  bf.Finalize();
  CHECK_EQUAL(12, bf.size());
  bf.Compute(0.0);
  Entity_ID_List face_list;
  mesh->get_set_entities(RIGHT, FACE, USED, &face_list);
  for (Entity_ID_List::iterator f = face_list.begin(); f != face_list.end(); ++f)
    CHECK_EQUAL(1.0, bf.find(*f)->second);
  mesh->get_set_entities(FRONT, FACE, USED, &face_list);
  for (Entity_ID_List::iterator f = face_list.begin(); f != face_list.end(); ++f)
    CHECK_EQUAL(2.0, bf.find(*f)->second);
  mesh->get_set_entities(BACK, FACE, USED, &face_list);
  for (Entity_ID_List::iterator f = face_list.begin(); f != face_list.end(); ++f)
    CHECK_EQUAL(3.0, bf.find(*f)->second);
}

TEST_FIXTURE(reference_mesh, values2)
{
  // Create the function f(t,x,y,z) = t * (x + 2y + 3z)
  std::vector<double> c(1,1.0);
  std::vector<int> p(1,1);
  std::auto_ptr<Function> f1(new PolynomialFunction(c, p));
  double g[3] = {1.0, 2.0, 3.0};
  std::vector<double> grad(g, g+3);
  std::auto_ptr<Function> f2(new LinearFunction(0.0, grad));
  Teuchos::RCP<CompositeFunction> f3 =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new SeparableFunction(f1,f2))));
  // Create the boundary function
  BoundaryFunction bf(mesh);
  std::vector<std::string> regions(2); regions[0] = RIGHT; regions[1] = BACK;
  bf.Define(regions, f3);
  // Check values at t=1
  bf.Compute(1.0);
  for (BoundaryFunction::Iterator i = bf.begin(); i != bf.end(); ++i) {
    AmanziGeometry::Point p = mesh->face_centroid(i->first);
    CHECK_EQUAL(p.x()+2*p.y()+3*p.z(), i->second);
  }
  // Check values at t=2
  bf.Compute(2.0);
  for (BoundaryFunction::Iterator i = bf.begin(); i != bf.end(); ++i) {
    AmanziGeometry::Point p = mesh->face_centroid(i->first);
    CHECK_EQUAL(2*(p.x()+2*p.y()+3*p.z()), i->second);
  }
}

TEST_FIXTURE(reference_mesh, bad_input)
{
  Teuchos::RCP<CompositeFunction> f =
    Teuchos::rcp(new CompositeFunction(Teuchos::rcp(new ConstantFunction(1.0))));
  BoundaryFunction bf(mesh);
  //bf.Define(INVALID, f); // no such face set
  CHECK_THROW(bf.Define(INVALID, f), Errors::Message);
  bf.Define(RIGHT, f);
  //bf.Define(RIGHT, f); // overlapping definition
  CHECK_THROW(bf.Define(RIGHT, f), Errors::Message);
}
