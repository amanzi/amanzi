# -*- mode: cmake -*-

#
# Functions for managing tests.
#

include(AddParallelTest)
include(CMakeParseArguments)
include(PrintVariable)



function(_ADD_TEST_LABELS test_name)

  get_test_property(${test_name} LABELS labels)
  if ("${labels}" STREQUAL "NOTFOUND")
    unset(labels)
  endif()

  list(APPEND labels "${ARGN}")
  list(REMOVE_DUPLICATES labels)
  set_tests_properties(${test_name} PROPERTIES LABELS "${labels}")

endfunction(_ADD_TEST_LABELS)

function (_REGISTER_TEST test_name test_exec nprocs labels mpi_args)

  foreach(nproc ${nprocs})
    if (${nproc} GREATER 1)
      ADD_PARALLEL_TEST(${test_name} ${test_exec} NPROCS ${nproc} MPI_EXEC_ARGS ${mpi_args})
      _add_test_labels(${test_name} "Parallel")
    else()
      add_test(${test_name} ${test_exec})
      _add_test_labels(${test_name} "Serial")
    endif()
  endforeach()

endfunction(_REGISTER_TEST)




function(ADD_UNIT_TEST test_name test_exec)

#
# Usage:
#
# ADD_UNIT_TEST(<test_name> <test_executable>
#                PARALLEL
#                [NPROCS procs1 ... ]
#                [MPI_EXEC_ARGS arg1 ... ]
#
# Option PARALLEL signifies that this is a parallel job. This is also
# implied by an NPROCS value > 1
#
# Optional NPROCS keyword starts a list of the number of processors to
# run the test on. Defaults to 1.
#
# Optional MPI_EXEC_ARGS keyword denotes extra arguments to give to
# mpi. It is ignored for serial tests.
#

  set(options "PARALLEL")
  set(oneValueArgs "")
  set(multiValueArgs NPROCS)

  cmake_parse_arguments(ADD_UNIT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(nprocs "${ADD_UNIT_TEST_NPROCS}" )
  if ("${nprocs}" STREQUAL "")
    set(nprocs 1)
  endif()

  set(labels ${ADD_UNIT_TEST_LABELS})
  set(mpi_args ${ADD_UNIT_TEST_MPI_EXEC_ARGS})

  _register_test("${test_name}" "${test_exec}" "${nprocs}" "${labels}" "${mpiargs}")
  _add_test_labels("${test_name}" "Unit")

endfunction(ADD_UNIT_TEST)




function(ADD_INTEGRATION_TEST test_name test_exec)

  set(options "PARALLEL")
  set(oneValueArgs "")
  set(multiValueArgs NPROCS FILES)

  cmake_parse_arguments(ADD_INTEGRATION_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(nprocs "${ADD_INTEGRATION_TEST_NPROCS}" )
  if ("${nprocs}" STREQUAL "")
    set(nprocs 1)
  endif()

  set(labels ${ADD_UNIT_TEST_LABELS})
  set(mpi_args ${ADD_UNIT_TEST_MPI_EXEC_ARGS})

  _register_test("${test_name}" "${test_exec}" "${nprocs}" "${labels}" "${mpiargs}")
  _add_test_labels(${test_name} "Integration")

  # TODO: Add dependencies on files

endfunction(ADD_INTEGRATION_TEST)


