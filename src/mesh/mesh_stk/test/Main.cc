#include <UnitTest++.h>
#include <TestReporterStdout.h>

#include <stk_util/parallel/Parallel.hpp>

stk::ParallelMachine parallel_machine;

int main(int argc, char *argv[])
{

    parallel_machine = stk::parallel_machine_init (&argc, &argv);

    int status = UnitTest::RunAllTests ();

    stk::parallel_machine_finalize ();

    return status;

}
