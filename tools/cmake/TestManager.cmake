# -*- mode: cmake -*-

#
# Functions for managing tests.
#

include(CMakeParseArguments)
include(PrintVariable)


function (_REGISTER_TEST test_name test_exec test_args nprocs is_parallel mpi_args)

  if ("${nprocs}" STREQUAL "")
    set(nprocs 1)
  endif()

      if ((${nprocs} GREATER 1) OR "${is_parallel}" OR "${TESTS_REQUIRE_MPIEXEC}")
      _add_parallel_test(${test_name} ${test_exec} "${test_args}" ${nprocs} "${mpi_args}")
      _add_test_labels(${test_name} "PARALLEL")
    else()
      add_test(${test_name} ${test_exec} ${test_args})
      _add_test_labels(${test_name} "SERIAL")
    ENDIF()

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
    _add_test_labels("${test_name}" "${match}")
  else()
    message(FATAL_ERROR, "No, or invalid test kind specified.")
  endif()

endfunction()


function(_build_mpiexec_command executable)

  # Check executable variable 
  if ( NOT executable )
    message(FATAL_ERROR "Must specify an executable in _build_mpiexec_command")
  endif()
  
  # Parse the remaning arguments
  set(optiotns     "")
  set(oneValueArgs "NPROCS;OUTPUT")
  set(multiValueArgs "EXEC_ARGS;MPIEXEC_ARGS")
  cmake_parse_arguments(PARSE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  # Need a return variable 
  print_variable(PARSE_OUTPUT)
  if ( NOT PARSE_OUTPUT )
    message(FATAL_ERROR "Must specify an output variable in _build_mpiexec_command")
  endif()

  if ( NOT PARSE_NPROCS )
    set(PARSE_NPROCS 1)
  endif()

  set(mpiexec_command "${MPI_EXEC} ${MPI_EXEC_NUMPROCS_FLAG} ${PARSE_NPROCS} ${MPI_EXEC_PREFLAGS}")
  foreach(arg ${PARSE_MPIEXEC_ARGS})
    set(mpiexec_command "${mpiexec_command} ${arg}")
  endforeach()

  set(exec_command_args)
  foreach (arg ${PARSE_EXEC_ARGS})
    set(exec_command_args "${exec_command_args} ${arg}")
  endforeach()

  set(${PARSE_OUTPUT} "${mpiexec_command} ${executable} ${MPI_EXEC_POSTFLAGS} ${exec_command_args}" PARENT_SCOPE)

endfunction(_build_mpiexec_command)



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

# Option PARALLEL signifies that this is a parallel job. This is also
# implied by an NPROCS value > 1
#
# Optional NPROCS keyword starts a list of the number of processors to
# run the test on. Defaults to 1.
#
# Optional MPI_EXEC_ARGS keyword denotes extra arguments to give to
# mpi. It is ignored for serial tests.


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

  # Labels, This is a CMake list type
  _add_test_kind_label(${test_name} ${AMANZI_TEST_KIND})
  get_test_property(${test_name} LABELS test_labels)
  if ( AMANZI_TEST_PARALLEL AND AMANZI_TEST_NPROCS )
    if ( ${AMANZI_TEST_NPROCS} GREATER 1 )
      list(APPEND test_labels PARALLEL)
    else()
      list(APPEND test_labels SERIAL)
    endif()  
  else()  
    list(APPEND test_labels SERIAL)
  endif()
  set_tests_properties(${test_name} PROPERTIES LABELS "${test_labels}")
  
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




