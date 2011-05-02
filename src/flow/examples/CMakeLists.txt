# -*- mode: cmake -*-

#
#  Amanzi
#   Flow examples
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(TestManager)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})
include_directories(${OUTPUT_SOURCE_DIR})
include_directories(${SOLVERS_SOURCE_DIR})
include_directories(${TIME_STEPPING_SOURCE_DIR})
include_directories(${MESH_FACTORY_SOURCE_DIR})

if ( ENABLE_MOAB_Mesh )

    # Binary: richards-nlk
    add_executable(richards-nlk richards-nlk.cc)
    target_link_libraries(richards-nlk solvers moab_mesh flow cell_geometry output)

    # Binary: richards-jfnk
    add_executable(richards-jfnk richards-jfnk.cc)
    target_link_libraries(richards-jfnk solvers moab_mesh flow cell_geometry output)

    # Binary: richards-bdf2
    add_executable(richards-bdf2 richards-bdf2.cc)
    target_link_libraries(richards-bdf2 solvers time_integration solvers mesh_factory stk_mesh flow cell_geometry output)

else()    

    message(SEND_ERROR "Can not build flow examples without enabling MOAB Mesh. "
                       "Enable Moab_Mesh by adding the following cmake options\n"
                        "-D ENABLE_MOAB_Mesh:BOOL=on -D MOAB_DIR:FILEPATH=<moab_install_Path")
endif()                    

