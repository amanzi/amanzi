# -*- mode: cmake -*-

#
# Amanzi
#
# Add a parallel test
#
#  Usage:
#
#  ADD_PARALLEL_TEST(<test_name> <test_exec> [arg1] [arg2] [arg3]
#                     [ NPROCS n ]
#                     [ MPI_EXEC_ARGS arg1 arg2 arg3 ] )
#  NPROCS == number of MPI ranks (number of processors)
#  MPI_EXEC_ARGS = Additional arguments for the MPI_EXEC binary
#
#  WARNING: This will not work with POE. Command order is not correct.


# CMake Modules
include(CMakeParseArguments)

macro(ADD_PARALLEL_TEST test_name test_exec mpi_args)
  set(_options  "")
  set(_oneValue "NPROCS")
  set(_multiValue "")
  cmake_parse_arguments(ADD_PARALLEL_TEST "${_options}" "${_oneValue}" "${_multiValue}" ${ARGN})

  # Default number of procs is 1
  set(_local_nprocs  1)
  if (ADD_PARALLEL_TEST_NPROCS )
    set(_local_nprocs ${ADD_PARALLEL_TEST_NPROCS})
  endif()

  # Do not allow mpi_exec to run more than the MAX_NPROCS
  if ( MPI_EXEC_MAX_NUMPROCS )
    if (${MPI_EXEC_MAX_NUMPROCS} LESS ${_local_nprocs}  )
      set(_local_nprocs ${MPI_EXEC_MAX_NUMPROCS})
      message(WARNING "Test ${test_name} requested too many procs. "
        "Will run ${test_name} with ${_local_nprocs} ranks")
    endif()
  endif()

  set(test_exec_args "${MPI_EXEC_ARGS_FLAG}" "${mpi_args}")

  # Build MPI_EXEC cpmmand
  # This will have to change if POE is used as MPI_EXEC!
  set(_mpi_cmd ${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} ${_local_nprocs})
  list(APPEND _mpi_cmd ${test_exec_args} ${test_exec})

  # Finally add the test
  add_test(${test_name} ${_mpi_cmd})

endmacro(ADD_PARALLEL_TEST)

