#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "Teuchos_GlobalMPISession.hpp"

#include "VerboseObject_objs.hh"
#include "Kokkos_Core.hpp"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  Kokkos::initialize(); 
  auto r =  UnitTest::RunAllTests ();
  Kokkos::finalize(); 
  return r; 
}

