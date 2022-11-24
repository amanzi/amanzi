#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include "VerboseObject_objs.hh"
#include "VerboseObject.hh"

int
main(int argc, char* argv[])
{
  int status = UnitTest::RunAllTests();
  return status;
}
