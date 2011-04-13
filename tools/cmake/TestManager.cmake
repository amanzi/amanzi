# -*- mode: cmake -*-

#
# Functions for managing tests.
#

include(AddParallelTest)
include(CMakeParseArguments)
include(PrintVariable)

function(ADD_UNIT_TEST test_name test_exec)
#
# Usage:
#
# REGISTER_UNIT_TEST(<test_name> <test_executable>
#                    [NPROCS <list of processor counts>]
#                    [LABELS <list of labels>]
#                    [MPI_EXEC_ARGS arg1 ... ]
#
# Optional NPROCS keyword starts a list of the number of processors to
# run the test on. Defaults to 1.
#
# Optional LABELS assign string labels to the test.
#
# Optional MPI_EXEC_ARGS keyword denotes extra arguments to give to mpi.
#

set(options "")
set(oneValueArgs "")
set(multiValueArgs NPROCS LABELS)

# Will define values:
#  REGISTER_UNIT_TEST_NPROCS
#  REGISTER_UNIT_TEST_LABELS
# And possibly...
# REGISTER_UNIT_TEST_UNPARSED_ARGUMENTS
cmake_parse_arguments(ADD_UNIT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

set(nprocs "${ADD_UNIT_TEST_NPROCS}" )
if (nprocs EQUAL "")
  set(nprocs 1)
endif()

set(labels ${ADD_UNIT_TEST_LABELS} unit)

message("Name: ${test_name}  Exec: ${test_exec}  NProcs: ${nprocs}  Labels: ${labels}")

foreach(nproc ${nprocs})
  if (${nproc} GREATER 1)
    ADD_PARALLEL_TEST(${test_name} ${test_exec} NPROCS ${nproc} MPI_EXEC_ARGS ${ADD_UNIT_TEST_MPI_EXEC_ARGS})
  else()
    add_test(${test_name} ${test_exec})
  endif()
endforeach()


endfunction(ADD_UNIT_TEST)

