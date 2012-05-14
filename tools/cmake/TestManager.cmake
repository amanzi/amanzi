# -*- mode: cmake -*-

#
# Functions for managing tests.
#

include(CMakeParseArguments)
include(PrintVariable)

function(_APPEND_TEST_LABEL test_name label)

  get_test_property(${test_name} LABELS current_labels)
  if (current_labels)
    set_tests_properties(${test_name} PROPERTIES LABELS "${current_labels};${label}")
  else()  
    set_tests_properties(${test_name} PROPERTIES LABELS "${label}")
  endif()

endfunction(_APPEND_TEST_LABEL)

function(_ADD_TEST_KIND_LABEL test_name kind_in)

  set(kind_prefixes UNIT INT REG)

  string(TOUPPER "${kind_in}" kind)

  foreach(kind_prefix ${kind_prefixes})
    string(REGEX MATCH "${kind_prefix}" match ${kind})
    if(match)
      break()
    endif()
  endforeach()

 if (match)
    _append_test_label(${test_name} ${match})
  else()
    message(FATAL_ERROR "Invalid test label ${kind_in} (Valid Labels:${kind_prefixes})")
  endif()

endfunction(_ADD_TEST_KIND_LABEL)


# Usage:
#
# ADD_AMANZI_TEST(<test_name> <test_executable>
#                  [arg1 ...]
#                  KIND [unit | int | reg]
#                  [SOURCE file1 file2  ...]
#                  [LINK_LIBS lib1 lib2 ...]
#                  [PARALLEL] [EXPECTED_FAIL]
#                  [NPROCS procs1 ... ]
#                  [MPI_EXEC_ARGS arg1 ... ])

#
# Arguments:
#  test_name: the name given to the resulting test in test reports
#  test_executable: The test executable which performs the test
#  arg1 ...: Additional arguments for the test executable
#
# Keyword KIND is required and should be one of unit, int or reg. 
#
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
# Optional SOURCE_FILES keyword that defines a list of source files
# required to build test_executable. An add_executable call will be made
# if this option is active.
#
# Optional LINK_LIBS keyword defines a list of link libraries or link options
# to link test_executable. An target_link_libraries call will be made if
# this option is active.

function(ADD_AMANZI_TEST test_name test_exec)

  # --- Initialize 

  # Check test_name
  if ( NOT test_name )
    message(FATAL_ERROR "Must define a test name.")
  endif()

  # Check test_exec 
  if ( NOT test_exec )
    message(FATAL_ERROR "Must specify test executable name")
  endif()

  # Parse through the remaining options
  set(options PARALLEL EXPECTED_FAIL)
  set(oneValueArgs KIND)
  set(multiValueArgs NPROCS SOURCE LINK_LIBS MPI_EXEC_ARGS)
  cmake_parse_arguments(AMANZI_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  set(test_args    "${AMANZI_TEST_UNPARSED_ARGUMENTS}")

  # --- Check options

  # Force each test to parallel run if mpiexec is required
  if(TESTS_REQUIRE_MPIEXEC)
    set(AMANZI_TEST_PARALLEL TRUE)
  endif()

  # Force each PARALLEL TRUE if NPROCS set 
  if(AMANZI_TEST_NPROCS AND ( "${AMANZI_TEST_NPROCS}" GREATER 1 ) )
    set(AMANZI_TEST_PARALLEL TRUE)
  endif()  

  # Default to nprocs=1 when running parallel
  if ( AMANZI_TEST_PARALLEL AND (NOT AMANZI_TEST_NPROCS) )
    set(AMANZI_TEST_NPROCS 1)
  endif() 

  # Test the value of number of procs value
  if(AMANZI_TEST_NPROCS)
    if(NOT ("${AMANZI_TEST_NPROCS}" GREATER 0) )
      message(FATAL_ERROR "${AMANZI_TEST_NPROCS} is an invalid NPROCS value.")
    endif()

    if(MPI_EXEC_MAX_NUMPROCS AND AMANZI_TEST_PARALLEL)
      if ( "${MPI_EXEC_MAX_NUMPROCS}" LESS "${AMANZI_TEST_NPROCS}")
        message(WARNING "Test ${test_name} request too many nprocs (${AMANZI_TEST_NPROCS}). "
                        "Will skip this test.")
        return()
      endif()
    endif()
  endif() 

  # Require a KIND value
  if ( NOT AMANZI_TEST_KIND )
    message(FATAL_ERROR "A test type has not been specified for ${test_name}.")
  endif()


  # --- Define the test executable

  # Add the source file definitions
  if(AMANZI_TEST_SOURCE)
    add_executable(${test_exec} ${AMANZI_TEST_SOURCE})
  endif()

  # Add link libraries
  if(AMANZI_TEST_LINK_LIBS)
    target_link_libraries(${test_exec} ${AMANZI_TEST_LINK_LIBS})
  endif()

  # --- Add test

  # Adjust the execuable name if NOT fullpath AND TESTS_REQUIRE_FULLPATH is set
  if ( TESTS_REQUIRE_FULLPATH )
    if ( NOT ("${test_exec}" MATCHES "^/") )
      set(_tmp      "${CMAKE_CURRENT_BINARY_DIR}/${test_exec}")
      set(test_exec "${_tmp}")
    endif()  
  endif()

  # Construct the test execution command
  set(add_test_exec)
  set(add_test_args)
  if (AMANZI_TEST_PARALLEL)

    if ( MPI_EXEC_GLOBAL_ARGS )
      separate_arguments(global_mpi_args UNIX_COMMAND "${MPI_EXEC_GLOBAL_ARGS}")
    endif() 

    set(add_test_exec ${MPI_EXEC})
    set(add_test_args
                      ${MPI_EXEC_NUMPROCS_FLAG}
                      ${AMANZI_TEST_NPROCS}
                      ${global_mpi_args}
                      ${AMANZI_TEST_MPI_EXEC_ARGS}
                      ${MPI_EXEC_PREFLAGS}
                      ${test_exec}
                      ${MPI_EXEC_POSTFLAGS}
                      ${test_args})
  else()
    set(add_test_exec ${test_exec})
    set(add_test_args ${test_args})
  endif()

  # Call add_test
  add_test(NAME ${test_name} COMMAND ${add_test_exec} ${add_test_args})

  # --- Add test properties

  # Labels
  _add_test_kind_label(${test_name} ${AMANZI_TEST_KIND})
  if ( AMANZI_TEST_PARALLEL AND AMANZI_TEST_NPROCS )
    if ( ${AMANZI_TEST_NPROCS} GREATER 1 )
      _append_test_label(${test_name} PARALLEL)
    else()
      _append_test_label(${test_name} SERIAL)
    endif()  
  else()  
    _append_test_label(${test_name} SERIAL)
  endif()
  
  # Remaining properties are single valued. Building 
  # test_properties as a list should get past the CMake parser.

  # Timeout
  if ( TESTS_TIMEOUT_THRESHOLD )
    list(APPEND test_properties TIMEOUT ${TESTS_TMIEOUT_THRESHOLD})
  endif()

  # CTest needs to know how procs this test needs
  if ( AMANZI_TEST_PARALLEL )
    list(APPEND test_properties PROCESSORS ${AMANZI_TEST_NPROCS})
  endif()

  # Set expected failure flag
  if ( AMANZI_TEST_EXPECTED_FAIL )
     list(APPEND test_properties WILL_FAIL TRUE)
  endif() 

  if ( test_properties )
    set_tests_properties(${test_name} PROPERTIES ${test_properties})
  endif()

endfunction(ADD_AMANZI_TEST)




