#  -*- mode: cmake -*-
#
# Build tool as a TPL: Spack
#  

# --- Define all the directories and common external project flags

message(STATUS ">>>>>>>> Build_Spack.cmake")

find_package(Git REQUIRED)

set(GIT ${GIT_EXECUTABLE})
set(SPACK_URL https://github.com/spack/spack.git)
set(XSDK_BRANCH master)

define_external_project_args(SPACK
    TARGET Spack
    BUILD_IN_SOURCE 0
    )

ExternalProject_Add(${SPACK_BUILD_TARGET}
    #DOWNLOAD_DIR ${TPL_INSTALL_PREFIX}
    DOWNLOAD_COMMAND cd ${TPL_INSTALL_PREFIX} COMMAND ${GIT} clone ${SPACK_URL} 

    CONFIGURE_COMMAND ls

    #BINARY_DIR ${TPL_INSTALL_PREFIX}/spack
    BUILD_COMMAND cd ${TPL_INSTALL_PREFIX}/spack COMMAND ${GIT} checkout ${XSDK_BRANCH} COMMAND ${GIT} pull
    BUILD_ALWAYS 1

    #INSTALL_DIR ${TPL_INSTALL_PREFIX}/spack
    INSTALL_COMMAND set(SPACK_BINARY ${TPL_INSTALL_PREFIX}/spack/bin/spack PARENT_SCOPE)
)

#ExternalProject_add_step(Spack CLONE_CD
#                         COMMAND cd ${TPL_INSTALL_PREFIX}
#            DEPENDERS download
#)
#ExternalProject_add_step(Spack CHECKOUT_CD
#                         COMMAND cd ${TPL_INSTALL_PREFIX}/spack
#            DEPENDERS build
#)
#ExternalProject_add_step(Spack PULL_CD
#                         COMMAND cd ${TPL_INSTALL_PREFIX}/spack
#            DEPENDERS install
#)

#set (SPACK_BINARY ${TPL_INSTALL_PREFIX}/spack/bin/spack)
message(STATUS ">>>>>>>> SPACK_BINARY: ${SPACK_BINARY}")

if ( FALSE )
# clone the repo
set(GIT_ARGS clone https://github.com/LLNL/spack.git)
execute_process(COMMAND git ${GIT_ARGS}
            WORKING_DIRECTORY ${TPL_INSTALL_PREFIX}
                RESULT_VARIABLE err_occurred 
                OUTPUT_VARIABLE SPACK_GIT_STATUS
                ERROR_VARIABLE err
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE
                )
if(err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
endif()

message(STATUS ">>>> JDM: SPACK_GIT_STATUS:      ${SPACK_GIT_STATUS}")
unset(GIT_ARGS)

if ( ENABLE_XSDK )
    # checkout the correct branch
    set(GIT_ARGS checkout barry/xsdk)
    execute_process(COMMAND git ${GIT_ARGS}
                WORKING_DIRECTORY ${TPL_INSTALL_PREFIX}/spack
                    RESULT_VARIABLE err_occurred 
                    OUTPUT_VARIABLE SPACK_GIT_STATUS
                    ERROR_VARIABLE err
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE
                    )
    if(err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
    endif()
    
    message(STATUS ">>>> JDM: SPACK_GIT_STATUS:      ${SPACK_GIT_STATUS}")
    unset(GIT_ARGS)
    
    # do a git pull
    set(GIT_ARGS pull)
    execute_process(COMMAND git ${GIT_ARGS}
                WORKING_DIRECTORY ${TPL_INSTALL_PREFIX}/spack
                    RESULT_VARIABLE err_occurred 
                    OUTPUT_VARIABLE SPACK_GIT_STATUS
                    ERROR_VARIABLE err
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE
                    )
    if(err_occurred)
    message(WARNING "Error executing git:\n ${cmd}\n${err}")
    set(cmd_output cmd_output-NOTFOUND)
    exit()
    endif()
    
    message(STATUS ">>>> JDM: SPACK_GIT_STATUS:      ${SPACK_GIT_STATUS}")
    unset(GIT_ARGS)
    
endif()
endif()

