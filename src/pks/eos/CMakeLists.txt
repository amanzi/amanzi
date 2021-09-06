# -*- mode: cmake -*-

#
#  Amanzi
#    Equations of state
#

#
# Define a project name
# After this command the following varaibles are defined
#   EOS_SOURCE_DIR
#   EOS_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(EOS)

# Amanzi include directories
include_directories(${ATK_SOURCE_DIR})

#
# Library: eos
#
set(eos_src_files 
    EOSDensityEvaluator.cc
    EOS_DensityIdealGas.cc
    EOS_DensityVaporInGas.cc
    EOS_DensityWater.cc
    EOS_DensityWaterFEHM.cc
    EOS_DensityWaterTabular.cc
    Henry_Solubility.cc
    IsobaricEOSEvaluator.cc
    LookupTable.cc
    MolarFractionGasEvaluator.cc
    ThermalConductivity_Water.cc
    SaturatedVaporPressure_Water.cc
    SaturatedVaporPressureFactory.cc
    ViscosityConstant.cc
    ViscosityWater.cc
    ViscosityWaterFEHM.cc
    ViscosityWaterTabular.cc
    ViscosityBaseFactory.cc
    ViscosityEvaluator.cc
   )

add_amanzi_library(eos
                   SOURCE ${eos_src_files} HEADERS ${eos_inc_files}
		   LINK_LIBS data_structures state atk ${Epetra_LIBRARIES})
if (APPLE AND BUILD_SHARED_LIBS)
  set_target_properties(eos PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()			  

#
# Install Targets
#
set(eos_inc_files 
    EOSDensityEvaluator.hh
    EOSDensityFactory.hh
    EOS_Density.hh
    EOS_DensityIdealGas.hh
    EOS_DensityVaporInGas.hh
    EOS_DensityWater.hh
    EOS_DensityWaterFEHM.hh
    EOS_DensityWaterTabular.hh
    eos_registration.hh
    Henry_Solubility.hh
    IsobaricEOSEvaluator.hh
    LookupTable.hh
    MolarFractionGasEvaluator.hh
    ThermalConductivity_Water.hh
    SaturatedVaporPressure.hh
    SaturatedVaporPressure_Water.hh
    SaturatedVaporPressureFactory.hh
    ViscosityBase.hh
    ViscosityBaseFactory.hh
    ViscosityConstant.hh
    ViscosityWater.hh
    ViscosityWaterFEHM.hh
    ViscosityWaterTabular.hh
    ViscosityEvaluator.hh)

# install headers into a subdirectory
install(FILES ${eos_inc_files} DESTINATION include/amanzi)

if (BUILD_TESTS) 
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${EOS_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${EOS_SOURCE_DIR} EQUAL ${EOS_BINARY_DIR}) )
        file(GLOB DataFiles "${EOS_SOURCE_DIR}/test/*.xml"
                            "${EOS_SOURCE_DIR}/test/*.eos")
        file(COPY ${DataFiles} DESTINATION ${EOS_BINARY_DIR}/test/)
    endif()

    set(amanzi_libs eos atk)

    # Test: Multiscale
    add_amanzi_test(eos_evaluators eos_evaluators
                    KIND unit
                    SOURCE test/Main.cc test/eos_evaluators.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

