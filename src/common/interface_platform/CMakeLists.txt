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
include_directories(${DBC_SOURCE_DIR})
include_directories(${ATK_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) include directories
include_directories(${Teuchos_INCLUDE_DIRS})
if ( CCSE_FOUND )
  include_directories(${CCSE_INCLUDE_DIRS})
endif()  

include_directories(${Boost_INCLUDE_DIRS})

if (ENABLE_Structured)
  add_definitions("-DENABLE_Structured")
  include_directories(${SAMR_SOURCE_DIR}) 
endif()

#
# interface library
#
# file(GLOB _hh_inc_files "*.hh")
set(interface_inc_files
  Simulator.hh
  InputConverter.hh
  InputConverterU.hh
  InputConverterU_Defs.hh
  InputParserIS_Defs.hh
  InputTranslator.hh)

add_amanzi_library(interface
                      SOURCE 
                           InputConverter.cc InputConverterU.cc
                           InputConverterU_State.cc InputConverterU_MPC.cc InputConverterU_PK.cc
                           InputConverterU_Flow.cc InputConverterU_Transport.cc
                           InputConverterU_Energy.cc InputConverterU_Chemistry.cc
                           InputConverterU_Solvers.cc InputConverterU_Mesh.cc InputConverterU_Output.cc
                           InputConverterS.cc
			   InputTranslator.cc
                      HEADERS ${interface_inc_files}
                      LINK_LIBS error_handling geometry ${Teuchos_LIBRARIES})

if (BUILD_TESTS) 
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${PLATFORM_INTERFACE_SOURCE_DIR})

    # Copy test directory files if an out of source build
    if (NOT (${PLATFORM_INTERFACE_SOURCE_DIR} EQUAL ${PLATFORM_INTERFACE_BINARY_DIR}) )
        execute_process(COMMAND ${CMAKE_COMMAND} -E 
          copy_directory ${PLATFORM_INTERFACE_SOURCE_DIR}/test ${PLATFORM_INTERFACE_BINARY_DIR}/test) 
    endif()

    set(amanzi_libs interface atk)

    # Test: transient Darcy flow 
    add_amanzi_test(converter_u_base converter_u_base
                    KIND int
                    SOURCE test/Main.cc test/converter_u_base.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES} ${Boost_LIBRARIES}
                                             ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})
 
    # Test: translate Amanzi-S input
    if (ENABLE_Structured)
      add_amanzi_test(converter_s converter_s
                      KIND int
                      SOURCE test/Main.cc test/converter_s.cc
                      LINK_LIBS ${amanzi_libs} ${CCSE_LIBRARIES} ${UnitTest_LIBRARIES} ${Boost_LIBRARIES}
                                               ${XERCES_LIBRARIES} ${XERCES_ICU_LIBRARIES})
    endif()
 
endif()


