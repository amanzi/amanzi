#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include <map>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "boundary-function.hh"
#include "constant-function.hh"
#include "polynomial-function.hh"
#include "linear-function.hh"
#include "separable-function.hh"
#include "function-factory.hh"
#include "errors.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

int main (int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests ();
}

TEST(empty)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  Teuchos::RCP<Mesh> mesh(factory(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  const std::map<int,double> &bfv = bf(0.0);
  CHECK_EQUAL(0, bfv.size());
}

TEST(basic)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  Teuchos::RCP<Mesh> mesh(factory(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  Teuchos::RCP<Function> f1(new ConstantFunction(1.0));
  Teuchos::RCP<Function> f2(new ConstantFunction(2.0));
  Teuchos::RCP<Function> f3(new ConstantFunction(3.0));
  int r[5] = {1, 2, 3, 4, 4};
  // Add a definition for a single side.
  {
    std::vector<int> reg(r, r+1);
    bf.Define(reg, f1);
    const std::map<int,double> &bfv = bf(0.0);
    CHECK_EQUAL(4,bfv.size());
  }
  // Add a definition for a couple other sides
  {
    std::vector<int> reg(r+1, r+3);
    bf.Define(reg, f1);
    const std::map<int,double> &bfv = bf(0.0);
    CHECK_EQUAL(12,bfv.size());
  }
  // Add a definition for yet more sides, but with duplicates.
  {
    std::vector<int> reg(r+3, r+5);
    bf.Define(reg, f1);
    const std::map<int,double> &bfv = bf(0.0);
    CHECK_EQUAL(16,bfv.size());
  }
}

TEST(values1)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  Teuchos::RCP<Mesh> mesh(factory(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  BoundaryFunction bf(mesh);
  Teuchos::RCP<Function> f1(new ConstantFunction(1.0));
  Teuchos::RCP<Function> f2(new ConstantFunction(2.0));
  Teuchos::RCP<Function> f3(new ConstantFunction(3.0));
  bf.Define(1, f1);
  bf.Define(2, f2);
  bf.Define(3, f3);
  const std::map<int,double> &bfv = bf(0.0);
  CHECK_EQUAL(12, bfv.size());
  Entity_ID_List face_list;
  mesh->get_set_entities(1, FACE, USED, &face_list);
  for (Entity_ID_List::iterator f = face_list.begin(); f != face_list.end(); ++f)
    CHECK_EQUAL(1.0, bfv.find(*f)->second);
  mesh->get_set_entities(2, FACE, USED, &face_list);
  for (Entity_ID_List::iterator f = face_list.begin(); f != face_list.end(); ++f)
    CHECK_EQUAL(2.0, bfv.find(*f)->second);
  mesh->get_set_entities(3, FACE, USED, &face_list);
  for (Entity_ID_List::iterator f = face_list.begin(); f != face_list.end(); ++f)
    CHECK_EQUAL(3.0, bfv.find(*f)->second);
}

TEST(values2)
{
  // Create a 2x2x2 mesh on [0,4]^3.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  Teuchos::RCP<Mesh> mesh(factory(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  // Create the function f(t,x,y,z) = t * (x + 2y + 3z)
  std::vector<double> c(1,1.0);
  std::vector<int> p(1,1);
  std::auto_ptr<Function> f1(new PolynomialFunction(c, p));
  double g[3] = {1.0, 2.0, 3.0};
  std::vector<double> grad(g, g+3);
  std::auto_ptr<Function> f2(new LinearFunction(0.0, grad));
  Teuchos::RCP<Function> f3(new SeparableFunction(f1,f2));
  // Create the boundary function
  BoundaryFunction bf(mesh);
  std::vector<int> regions(2); regions[0] = 1; regions[1] = 3;
  // Check values at t=1
  for (BoundaryFunction::Iterator i = bf.begin(); i != bf.end(); ++i) {
    AmanziGeometry::Point p = mesh->face_centroid(i->first);
    CHECK_EQUAL(p.x()+2*p.y()+3*p.z(), i->second);
  }
  // Check values at t=2
  for (BoundaryFunction::Iterator i = bf.begin(); i != bf.end(); ++i) {
    AmanziGeometry::Point p = mesh->face_centroid(i->first);
    CHECK_EQUAL(2*(p.x()+2*p.y()+3*p.z()), i->second);
  }
}

TEST(bad_input)
{
  // Create a 2x2x2 mesh on [0,4]^3.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  MeshFactory factory(comm);
  Teuchos::RCP<Mesh> mesh(factory(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2));
  Teuchos::RCP<Function> f(new ConstantFunction(1.0));
  BoundaryFunction bf(mesh);
  //bf.Define(99, f); // no such face set
  CHECK_THROW(bf.Define(99, f), Errors::Message);
  bf.Define(1, f);
  //bf.Define(1, f); // overlapping definition
  CHECK_THROW(bf.Define(1, f), Errors::Message);
}
