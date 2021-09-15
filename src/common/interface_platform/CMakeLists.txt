# -*- mode: cmake -*-

#
#  Amanzi
#    Interface Platform
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   PLATFORM_INTERFACE_SOURCE_DIR
#   PLATFORM_INTERFACE_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(PLATFORM_INTERFACE)

# Amanzi include directories
include_directories(${AMANZI_BINARY_DIR}) # required to pick up amanzi_version.hh
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
if (CCSE_FOUND)
  include_directories(${CCSE_INCLUDE_DIRS})
endif()  

include_directories(${Boost_INCLUDE_DIRS})

#
# interface library
#
# file(GLOB _hh_inc_files "*.hh")
set(interface_inc_files
  Simulator.hh
  InputConverter.hh
  InputConverterU.hh
  InputConverterU_Defs.hh)

set(interface_src_files InputConverter.cc InputConverterU.cc
                        InputConverterU_State.cc InputConverterU_PK.cc
                        InputConverterU_MPC.cc InputConverterU_MPC_deprecated.cc
                        InputConverterU_Flow.cc InputConverterU_Transport.cc
                        InputConverterU_Energy.cc
                        InputConverterU_Chemistry.cc InputConverter_ThermoDatabase.cc
                        InputConverterU_ShallowWater.cc
                        InputConverterU_Multiphase.cc
                        InputConverterU_Solvers.cc InputConverterU_Mesh.cc
                        InputConverterU_Output.cc)

add_amanzi_library(interface 
                   SOURCE ${interface_src_files}
                   HEADERS ${interface_inc_files}
                   LINK_LIBS error_handling geometry atk
		             ${Teuchos_LIBRARIES} ${Boost_LIBRARIES}
                             ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})

if (ENABLE_Structured)
  add_definitions("-DENABLE_Structured")
  include_directories(../../structured_grid/mpc_pk) 
  link_directories(${CCSE_LIBRARY_DIR})

  add_amanzi_library(interface_s
                     SOURCE InputConverterS.cc
                     HEADERS ${interface_inc_files}
		     LINK_LIBS error_handling geometry interface
		     ${Teuchos_LIBRARIES} ${Boost_LIBRARIES} ${CCSE_LIBRARIES}
                               ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})
endif()

if (BUILD_TESTS) 
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${PLATFORM_INTERFACE_SOURCE_DIR} EQUAL ${PLATFORM_INTERFACE_BINARY_DIR}) )
      file(GLOB DataFiles "${PLATFORM_INTERFACE_SOURCE_DIR}/test/*.xml"
                          "${PLATFORM_INTERFACE_SOURCE_DIR}/test/*.dat"
                          "${PLATFORM_INTERFACE_SOURCE_DIR}/test/*.dump")
      file(COPY ${DataFiles} DESTINATION ${PLATFORM_INTERFACE_BINARY_DIR}/test/)
    endif()

    # Test: transient Darcy flow 
    add_amanzi_test(converter_u_base converter_u_base
                    KIND int
                    SOURCE test/Main.cc test/converter_u_base.cc 
                    LINK_LIBS interface atk 
                              ${UnitTest_LIBRARIES} ${Boost_LIBRARIES}
                              ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})
 
    # Test: translate Amanzi-S input
    if (ENABLE_Structured)
      if (BUILD_SHARED_LIBS)
        set(closing_dso test/converter_s_dso.f90)
      endif()

      add_amanzi_test(converter_s converter_s
                      KIND int
                      SOURCE test/Main.cc test/converter_s.cc ${closing_dso}
                      LINK_LIBS interface interface_s atk 
                                ${UnitTest_LIBRARIES} ${Boost_LIBRARIES}
                                ${CCSE_LIBRARIES} ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})
    endif()
endif()


