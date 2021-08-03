# ############################################################################ #
#
# DefineCompilerVersion
#  
# ############################################################################ #

include(PrintVariable)

function(DEFINE_COMPILER_VERSION)

   # Define languages
   set(_languages C CXX Fortran)
   foreach ( _lang IN LISTS _languages)

       # Define the full path name of the compiler and the id
       # Default value for each language is NOTFOUND
       set(_lang_compiler ${CMAKE_${_lang}_COMPILER})
       set(_compiler_id   ${CMAKE_${_lang}_COMPILER_ID})
       set(_version       CMAKE_${_lang}_COMPILER_VERSION-NOTFOUND)

       message(STATUS "Identifying ${_lang_compiler} (${_compiler_id}) version")

       # Check only if the ID is set
       if ( _compiler_id )
         string(TOUPPER ${_compiler_id} _compiler_id_uc)
         set(_version_cmd_opt)
         set(_regexp_pattern)

         # For now, I assume the version option is the same for all languages in
         # a compiler group. This may need tweaking in the future.
         if ( ${_compiler_id_uc} STREQUAL "GNU" )
           set(_version_cmd_opt "--version")
           set(_regexp_pattern ".*\\(.*\\)[ ]+([0-9]+\\.[0-9]+\\.[0-9]+).*")	
         elseif(${_compiler_id_uc} STREQUAL "PGI" )
           set(_version_cmd_opt "-V")
           set(_regexp_pattern ".*pgcc[ ]+([0-9]+\\.[0-9]+-[0-9]+).*")
         elseif (${_compiler_id_uc} STREQUAL "INTEL" )
           set(_version_cmd_opt "-v")
	   set(_regexp_pattern ".*,[ ]+Version ([0-9]+\\.[0-9]+\\.*[0-9]*).*Build.*")
         elseif (${_compiler_id_uc} STREQUAL "PATHSCALE" )
           set(_version_cmd_opt "--version")
           set(_regexp_pattern ".*Version ([0-9]+\\.[0-9]+\\.[0-9]+).*")
         elseif (${_compiler_id_uc} STREQUAL "CRAY" )
           set(_version_cmd_opt "-V")
           set(_regexp_pattern ".*Cray[ ]+C[ ]+:[ ]+Version ([0-9]+\\.[0-9]+\\.[0-9]+).*")
         elseif (${_compiler_id_uc} STREQUAL "CLANG" )
           set(_version_cmd_opt "-v")
           set(_regexp_pattern ".*\(clang-\)[ ]+([0-9]+\\.[0-9]+\\.[0-9]+).*")
         elseif (${_compiler_id_uc} STREQUAL "APPLECLANG" )
           set(_version_cmd_opt "-v")
           set(_regexp_pattern ".*\(clang-\)[ ]+([0-9]+\\.[0-9]+\\.[0-9]+).*")
         else()
           message(WARNING "Unknown compiler ID type ${_compiler_id_uc}")
         endif()

         # Execute the command if the option was set
         if (DEFINED _version_cmd_opt)
           execute_process(COMMAND ${_lang_compiler} ${_version_cmd_opt}
                           RESULT_VARIABLE result
                           OUTPUT_VARIABLE output
                           ERROR_VARIABLE  output 
                           OUTPUT_STRIP_TRAILING_WHITESPACE
                           ERROR_STRIP_TRAILING_WHITESPACE)

           # result > 0 indicates an error return 
           if ( result ) 
             message(SEND_ERROR "${_lang_compiler} ${_version_cmd_opt} failed."
                                "Ouptut:\n${output}")
           else()
             # string(REGEX REPLACE "([^\n]+).*" "\\1" first_line "${output}")
	     # string(REGEX REPLACE "${_regexp_pattern}" "\\1" _version "${first_line}")
             string(REGEX REPLACE "${_regexp_pattern}" "\\1" _version "${output}")
             # A fix for the PGI case because they use a '-' to separate the patch number. Annoying.
             if ( ${_compiler_id_uc} STREQUAL "PGI")
               set(_tmp ${_version})
               string(REGEX REPLACE "-" "." _version ${_tmp})
             endif()
	     string(STRIP "${_version}" _version)
           endif()  

         endif()

       else() 
         message(SEND_ERROR "The ${_lang} compiler ID is not defined")
       endif()

       # Message sent to indicate version definition
       message(STATUS "Identifying ${_lang_compiler} (${_compiler_id}) version -- ${_version}")

       # Push the definition back to the calling routine
	   set(CMAKE_${_lang}_COMPILER_VERSION ${_version} PARENT_SCOPE)

   endforeach()    


endfunction(DEFINE_COMPILER_VERSION) 

