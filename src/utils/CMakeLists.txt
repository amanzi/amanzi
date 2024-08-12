# -*- mode: cmake -*-
#
#  Amanzi Tool Kit (ATK)
#    General utilities
#    This directory has become a dumping ground for things that do not
#    have a well-defined home. Will go away or get renamed. For now we
#    need to keep it around.
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)
include(LibraryManager)
include(BinaryManager)

#
# Define a project name
# After this command the following varaibles are defined
#   ATK_SOURCE_DIR
#   ATK_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(ATK)

include_directories(${EXPRTK_INLUCDE_DIR})
include_directories(${Teuchos_INCLUDE_DIRS})
include_directories(${Epetra_INCLUDE_DIRS})

include_directories(${DBC_SOURCE_DIR})

message(STATUS "ETC: -- linking NetCDF_C_LIBRARIES = ${NetCDF_C_LIBRARIES}")

set(tpl_link_libs
  ${Epetra_LIBRARIES}
  ${Teuchos_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${NetCDF_C_LIBRARIES}
)

add_amanzi_library(atk
                   SOURCE AmanziComm.cc
                          ExprTK.cc
                          Formatter.cc
                          ReaderHDF5.cc
                          ReaderNetCDF.cc
                          IOEvent.cc
                          Key.cc
                          LeastSquare.cc 
                          Spline.cc
                          SplinedCurve.cc
                          TimeStepManager.cc
                          Units.cc
                          VerboseObject.cc
                          XMLParameterListWriter.cc 
                   HEADERS Abbreviations.hh
                           AmanziComm.hh
                           AmanziTypes.hh
                           Brent.hh
                           ExprTK.hh
                           Factory.hh
                           FIFO_Set.hh
                           Formatter.hh
                           Reader.hh
                           ReaderHDF5.hh
                           ReaderNetCDF.hh
                           IOEvent.hh
                           Iterators.hh
                           Key.hh
                           MeshView.hh
                           LeastSquare.hh
                           ObservationData.hh
                           Spline.hh
                           SplinedCurve.hh
                           StringExt.hh
                           Tag.hh
                           TimeStepManager.hh
                           Units.hh
                           Utils.hh
                           VerboseObject.hh
                           VerboseObject_objs.hh
                           ViewUtils.hh
                           XMLParameterListWriter.hh
                           StringReducer.hh
                    LINK_LIBS error_handling ${tpl_link_libs})

if (APPLE AND BUILD_SHARED_LIBS)
  set_target_properties(atk PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()			  

add_amanzi_executable(xml_to_yaml
  SOURCE XmlToYaml.cc
  LINK_LIBS ${tpl_link_libs}
  OUTPUT_NAME xml_to_yaml
  OUTPUT_DIRECTORY ${ATS_BINARY_DIR})

add_amanzi_executable(yaml_to_xml
  SOURCE YamlToXml.cc
  LINK_LIBS ${tpl_link_libs}
  OUTPUT_NAME yaml_to_xml
  OUTPUT_DIRECTORY ${ATS_BINARY_DIR})


if (BUILD_TESTS)
    # Add UnitTest include directoy
    include_directories(${UnitTest_INCLUDE_DIRS})

    # Copy test directory files if an out of source build
    if (NOT (${ATK_SOURCE_DIR} EQUAL ${ATK_BINARY_DIR}))
      file(GLOB DataFiles "${ATK_SOURCE_DIR}/test/*")
      file(COPY ${DataFiles} DESTINATION ${ATK_BINARY_DIR}/test/)
    endif()

    # Add the utils directory to the include paths
    include_directories(${ATK_SOURCE_DIR})

    add_amanzi_test(utils_keys utils_keys
                    KIND unit
                    SOURCE test/Main.cc test/utils_keys.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})
    
    add_amanzi_test(utils_time_step_manager utils_time_step_manager
                    KIND unit
                    SOURCE test/Main.cc test/utils_time_step_manager.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_verbosity utils_verbosity
                    KIND unit
                    SOURCE test/Main.cc test/utils_verboseobject.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_formatter utils_formatter
                    KIND unit
                    SOURCE test/Main.cc test/utils_formatter.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_units utils_units
                    KIND unit
                    SOURCE test/Main.cc test/utils_units.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_spline utils_spline
                    KIND unit
                    SOURCE test/Main.cc test/utils_spline.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_ioevent utils_ioevent
                    KIND unit
                    SOURCE test/Main.cc test/utils_ioevent.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_exprtk utils_exprtk
                    KIND unit
                    SOURCE test/Main.cc test/utils_exprtk.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_reader utils_reader
                    KIND unit
                    SOURCE test/Main.cc test/utils_reader.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

    add_amanzi_test(utils_reducer utils_reducer
                    KIND unit
                    SOURCE test/Main.cc test/utils_reducer.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})
    add_amanzi_test(utils_reducer_np2 utils_reducer NPROCS 2 KIND unit)
    add_amanzi_test(utils_reducer_np5 utils_reducer NPROCS 5 KIND unit)

    add_amanzi_test(utils_brent utils_brent
                    KIND unit
                    SOURCE test/Main.cc test/utils_brent.cc
                    LINK_LIBS atk ${UnitTest_LIBRARIES})

endif()

