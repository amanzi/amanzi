# -*- mode: cmake -*-

#
# Functions for managing tests.
#

include(AddParallelTest)
include(CMakeParseArguments)
include(PrintVariable)


function (_REGISTER_TEST test_name test_exec test_args nprocs is_parallel mpi_args)

  if ("${nprocs}" STREQUAL "")
    set(nprocs 1)
  endif()

  foreach(nproc ${nprocs})
    if ((${nproc} GREATER 1) OR "${is_parallel}")
      _add_parallel_test(${test_name} ${test_exec} "${test_args}" ${nproc} "${mpi_args}")
      _add_test_labels(${test_name} "Parallel")
    else()
      add_test(${test_name} ${test_exec} ${test_args})
      _add_test_labels(${test_name} "Serial")
    endif()
  endforeach()

endfunction(_REGISTER_TEST)



function(_ADD_PARALLEL_TEST test_name test_exec test_args nproc mpi_args)

  set(_options  "")
  set(_oneValue "")
  set(_multiValue "MPI_EXEC_ARGS")
  cmake_parse_arguments(ADD_PARALLEL_TEST "${_options}" "${_oneValue}" "${_multiValue}" ${ARGN}) 

  # Do not allow mpi_exec to run more than the MAX_NPROCS
  if ( MPI_EXEC_MAX_NUMPROCS )
    if (${MPI_EXEC_MAX_NUMPROCS} LESS ${nproc}  )
      set(nproc ${MPI_EXEC_MAX_NUMPROCS})
      message(WARNING "Test ${test_name} requested too many procs. "
        "Will run ${test_name} with ${nproc} ranks")
    endif()
  endif()

  # Build MPI_EXEC cpmmand
  # This will have to change if POE is used as MPI_EXEC!
  set(_mpi_cmd "${MPI_EXEC}" "${MPI_EXEC_NUMPROCS_FLAG}" "${nproc}" "${MPI_EXEC_ARGS_FLAG}" "${mpi_args}" "${test_exec}" "${test_args}")

  # Register the test
  add_test(${test_name} ${_mpi_cmd})

endfunction(_ADD_PARALLEL_TEST)




function(_ADD_TEST_LABELS test_name)

  get_test_property(${test_name} LABELS labels)
  if ("${labels}" STREQUAL "NOTFOUND")
    unset(labels)
  endif()

  list(APPEND labels "${ARGN}")
  list(REMOVE_DUPLICATES labels)
  set_tests_properties(${test_name} PROPERTIES LABELS "${labels}")

endfunction(_ADD_TEST_LABELS)








#
# Usage:
#
# ADD_UNIT_TEST(<test_name> <test_executable>
#                [arg1 ...]
#                PARALLEL
#                [NPROCS procs1 ... ]
#                [MPI_EXEC_ARGS arg1 ... ]
#
# Arguments:
#  test_name: the name given to the resulting test in test reports
#  test_executable: The test executable which performs the test
#  arg1 ...: Additional arguments for the test executable
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

function(ADD_UNIT_TEST test_name test_exec)

  set(options "PARALLEL")
  set(oneValueArgs "")
  set(multiValueArgs NPROCS)

  cmake_parse_arguments(ADD_UNIT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(test_args    "${ADD_UNIT_TEST_UNPARSED_ARGUMENTS}")
  set(is_parallel  "${ADD_UNIT_TEST_PARALLEL}")
  set(mpi_args     "${ADD_UNIT_TEST_MPI_EXEC_ARGS}")
  set(nprocs       "${ADD_UNIT_TEST_NPROCS}")

  separate_arguments(global_mpi_args UNIX_COMMAND "${MPI_EXEC_ARGS}")
  list(APPEND mpi_args ${global_mpi_args})

  _register_test("${test_name}" "${test_exec}" "${test_args}" "${nprocs}" "${is_parallel}" "${mpi_args}")
  _add_test_labels("${test_name}" "Unit")

endfunction(ADD_UNIT_TEST)


function(ADD_INTEGRATION_TEST test_name test_exec)

  set(options "PARALLEL")
  set(oneValueArgs "")
  set(multiValueArgs NPROCS FILES)

  cmake_parse_arguments(ADD_INTEGRATION_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  set(test_args    "${ADD_INTEGRATION_TEST_UNPARSED_ARGUMENTS}")
  set(is_parallel  "${ADD_INTEGRATION_TEST_PARALLEL}")
  set(mpi_args     "${ADD_INTEGRATION_TEST_MPI_EXEC_ARGS}")
  set(nprocs       "${ADD_INTEGRATION_TEST_NPROCS}")

  _register_test("${test_name}" "${test_exec}" "${test_args}" "${nprocs}" "${is_parallel}" "${mpi_args}")
  _add_test_labels(${test_name} "Integration")

  # TODO: Add dependencies on files

endfunction(ADD_INTEGRATION_TEST)


