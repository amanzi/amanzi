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


