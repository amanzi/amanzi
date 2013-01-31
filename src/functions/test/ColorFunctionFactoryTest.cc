#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_MpiComm.h"

#include <string>

#include "ColorFunction.hh"
#include "ColorFunctionFactory.hh"
#include "errors.hh"

using namespace Amanzi;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests();
}

TEST(Grid1D_cell)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  ColorFunctionFactory factory;
  std::string infile("test/ColorFunc1CellTest.txt");
  ColorFunction *f = factory.Create(infile, comm);
  double x;
  x = 0.0;
  CHECK_EQUAL(10, (*f)(&x));
  x = 1.0;
  CHECK_EQUAL(20, (*f)(&x));
}

TEST(Grid1D_node)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  ColorFunctionFactory factory;
  std::string infile("test/ColorFunc1NodeTest.txt");
  ColorFunction *f = factory.Create(infile, comm);
  double x;
  x = 0.0;
  CHECK_EQUAL(10, (*f)(&x));
  x = 1.0;
  CHECK_EQUAL(20, (*f)(&x));
}

TEST(Grid2D_cell)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  ColorFunctionFactory factory;
  std::string infile("test/ColorFunc2CellTest.txt");
  ColorFunction *f = factory.Create(infile, comm);
  double x1[2] = {2.0, 1.0};
  CHECK_EQUAL(30, (*f)(x1));
  double x2[2] = {0.0, 2.0};
  CHECK_EQUAL(40, (*f)(x2));
  double x3[2] = {0.0, 1.0};
  CHECK_EQUAL(10, (*f)(x3));
  double x4[2] = {2.0, 2.0};
  CHECK_EQUAL(60, (*f)(x4));
}

TEST(Grid2D_node)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  ColorFunctionFactory factory;
  std::string infile("test/ColorFunc2NodeTest.txt");
  ColorFunction *f = factory.Create(infile, comm);
  double x1[2] = {2.0, 1.0};
  CHECK_EQUAL(30, (*f)(x1));
  double x2[2] = {0.0, 2.0};
  CHECK_EQUAL(40, (*f)(x2));
  double x3[2] = {0.0, 1.0};
  CHECK_EQUAL(10, (*f)(x3));
  double x4[2] = {2.0, 2.0};
  CHECK_EQUAL(60, (*f)(x4));
}

TEST(Grid3D_cell)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  ColorFunctionFactory factory;
  std::string infile("test/ColorFunc3CellTest.txt");
  ColorFunction *f = factory.Create(infile, comm);
  // check some random locations
  double x1[3] = {2.0, 1.0, 3.0};
  CHECK_EQUAL(213, (*f)(x1));
  double x2[3] = {2.0, 2.0, 0.0};
  CHECK_EQUAL(220, (*f)(x2));
  double x3[3] = {3.0, 1.0, 1.0};
  CHECK_EQUAL(311, (*f)(x3));
  double x4[3] = {3.0, 3.0, 3.0};
  CHECK_EQUAL(333, (*f)(x4));
}

TEST(Grid3D_node)
{
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  ColorFunctionFactory factory;
  std::string infile("test/ColorFunc3NodeTest.txt");
  ColorFunction *f = factory.Create(infile, comm);
  // check some random locations
  double x1[3] = {2.0, 1.0, 3.0};
  CHECK_EQUAL(213, (*f)(x1));
  double x2[3] = {2.0, 2.0, 0.0};
  CHECK_EQUAL(220, (*f)(x2));
  double x3[3] = {3.0, 1.0, 1.0};
  CHECK_EQUAL(311, (*f)(x3));
  double x4[3] = {3.0, 3.0, 3.0};
  CHECK_EQUAL(333, (*f)(x4));
}
