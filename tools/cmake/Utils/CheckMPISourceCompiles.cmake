#
# CHECK_MPI_SOURCE_COMPILES
#
#  A function that verifies that the CMAKE_C_COMPILER, 
#  CMAKE_CXX_COMPILER and CMAKE_Fortran_COMPILER can compile and link 
#  a MPI test file. Since this function is a simple wrapper
#  call to check_c{xx}_source_compiles macros, the following
#  variables can be set before the this function to alter
#  how the mocros behave:
#         CMAKE_REQUIRED_FLAGS = string of compile command line flags
#         CMAKE_REQUIRED_DEFINITIONS = list of macros to define (-DFOO=bar)
#         CMAKE_REQUIRED_INCLUDES = list of include directories
#         CMAKE_REQUIRED_LIBRARIES = list of libraries to link
#
function(CHECK_MPI_SOURCE_COMPILES FLAG)

  include(CheckCSourceCompiles)

  message(STATUS "Checking whether C compiler can compile MPI program")

  set(_mpi_c_source "
      #include <stdio.h>
      #include <mpi.h>
      void main(int argc, char **argv) 
      {
        MPI_Init(&argc,&argv);
        puts(__LINE__);
        MPI_Finalize();
      }")

    # CHECK_C_SOURCE_COMPILES is a MACRO not a function
    # check $VAR for compile status.
    check_c_source_compiles("${_mpi_c_source}" MPI_C)
    if(${MPI_C})
      message(STATUS "Checking whether C compiler can compile MPI program - yes")
    else()
      message(STATUS "Checking whether C compiler can compile MPI program - no")
    endif()   
        


    include(CheckCXXSourceCompiles)
    message(STATUS "Checking whether C++ compiler can compile MPI program")
    set(_mpi_cxx_source "
        #include <iostream>
        #include \"mpi.h\"
        int main(int argc, char *argv[])
        {
          MPI::Init(argc,argv);
          MPI::Finalize();
          return 0;
        }")
    
    # CHECK_C_SOURCE_COMPILES is a MACRO not a function
    # check $VAR for compile status.
    check_cxx_source_compiles("${_mpi_cxx_source}" MPI_CXX)

    if(${MPI_CXX})
      message(STATUS "Checking whether C++ compiler can compile MPI program - yes")
    else()
      message(STATUS "Checking whether C++ compiler can compile MPI program - no")
    endif() 

    # Check the Fortran compiler if enabled
    if ( CMAKE_Fortran_COMPILER ) 
      message(STATUS "Checking whether Fortran compiler can compile MPI program")
      set(_mpi_fortran_source "
          PROGRAM mpi_test
          INTEGER ierr
          call MPI_Init(ierr)
          call MPI_Finalize(ierr)
          STOP
          END PROGRAM")

      # As of CMAKE 2.8.6 no CheckFortranSourceCompiles macro
      set(test_file "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_mpi.f") 
      file(WRITE "${test_file}" "${_mpi_fortran_source}\n")

      message(STATUS "Performing Test MPI_Fortran")
      try_compile(MPI_Fortran
                  ${CMAKE_BINARY_DIR}
                  ${test_file}
                  OUTPUT_VARIABLE result)
      if(${MPI_Fortran})
        message(STATUS "Performing Test MPI_Fortran - Success")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
                    "Performing Fortran SOURCE FILE Test MPI_Fortran succeeded with the following output:\n"
                    "${result}\n"
                    "Source file was:\n${_mpi_fortran_source}\n")
        message(STATUS "Checking whether Fortran compiler can compile MPI program - yes")
      else()
        message(STATUS "Performing Test MPI_Fortran - Failed")
        set(MPI_Fortran "" CACHE INTERNAL "Test MPI_Fortran")
        file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
                    "Performing Fortran SOURCE FILE Test MPI_Fortran failed with the following output:\n"
                    "${result}\n"
                    "Source file was:\n${_mpi_fortran_source}\n")
        message(STATUS "Checking whether Fortran compiler can compile MPI program - no")
       endif()            

     endif()

     # Now push up the flag to the parent namespace
     if ( CMAKE_Fortran_COMPILER )
       if ( "${MPI_C}" AND "${MPI_CXX}" AND "${MPI_Fortran}" )
         set(${FLAG} TRUE PARENT_SCOPE)
       else()
         set(${FLAG} FALSE PARENT_SCOPE)
       endif()
     else()
       if ( "${MPI_C}" AND "${MPI_CXX}" )
         set(${FLAG} TRUE PARENT_SCOPE)
       else()
         set(${FLAG} FALSE PARENT_SCOPE)
       endif()
     endif()  
         
endfunction(CHECK_MPI_SOURCE_COMPILES)
