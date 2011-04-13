# -*- mode: cmake -*-

#
# Functions for managing tests.
#

include(CMakeParseArguments)
include(PrintVariable)

function(REGISTER_UNIT_TEST test_name test_exec)
#
# Usage:
#
# REGISTER_UNIT_TEST(<test_name> <test_executable>
#                    [NPROCS <list of processor counts>]
#                    [LABELS <list of labels>]
#
# Optional NPROCS keyword starts a list of the number of processors to
# run the test on. Defaults to 1.
#
# Optional LABELS assign string labels to the test.

set(options "")
set(oneValueArgs "")
set(multiValueArgs NPROCS LABELS)

# Will define values:
#  REGISTER_UNIT_TEST_NPROCS
#  REGISTER_UNIT_TEST_LABELS
# And possibly...
# REGISTER_UNIT_TEST_UNPARSED_ARGUMENTS
cmake_parse_arguments(REGISTER_UNIT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

set(nprocs "${REGISTER_UNIT_TEST_NPROCS}" )
if (nprocs EQUAL "")
  set(nprocs 1)
endif()

set(labels ${REGISTER_UNIT_TEST_LABELS} unit)

message("Name: ${test_name}  Exec: ${test_exec}  NProcs: ${nprocs}  Labels: ${labels}")

endfunction(REGISTER_UNIT_TEST)
