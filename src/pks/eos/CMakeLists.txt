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
    EOSEvaluator.cc
    EOSFactory.cc 
    EOS_Constant.cc
    EOS_IdealGas.cc
    EOS_VaporInGas.cc
    EOS_Water.cc
    EOS_WaterFEHM.cc
    Henry_Solubility.cc
    IsobaricEOSEvaluator.cc
    MolarFractionGasEvaluator.cc
    ThermalConductivity_Water.cc
    SaturatedVaporPressure_Water.cc
    SaturatedVaporPressureFactory.cc
    Viscosity_Water.cc
    Viscosity_Constant.cc
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
    EOS.hh
    EOSEvaluator.hh
    EOSFactory.hh
    EOS_Constant.hh
    EOS_ConstantMolarMass.hh
    EOS_IdealGas.hh
    EOS_VaporInGas.hh
    EOS_Water.hh
    EOS_WaterFEHM.hh
    eos_registration.hh
    Henry_Solubility.hh
    IsobaricEOSEvaluator.hh
    MolarFractionGasEvaluator.hh
    ThermalConductivity_Water.hh
    SaturatedVaporPressure.hh
    SaturatedVaporPressure_Water.hh
    SaturatedVaporPressureFactory.hh
    Viscosity_Base.hh
    Viscosity_Constant.hh
    Viscosity_Water.hh
    ViscosityBaseFactory.hh
    ViscosityEvaluator.hh)

# install headers into a subdirectory
install(FILES ${eos_inc_files} DESTINATION include/amanzi)

if (BUILD_TESTS) 
    include_directories(${UnitTest_INCLUDE_DIRS})
    include_directories(${EOS_SOURCE_DIR})

    # Copy test directory files if any out of source build
    if (NOT (${EOS_SOURCE_DIR} EQUAL ${EOS_BINARY_DIR}) )
        file(GLOB DataFiles "${EOS_SOURCE_DIR}/test/*.xml")
        file(COPY ${DataFiles} DESTINATION ${EOS_BINARY_DIR}/test/)
    endif()

    set(amanzi_libs eos atk)

    # Test: Multiscale
    add_amanzi_test(eos_evaluators eos_evaluators
                    KIND unit
                    SOURCE test/Main.cc test/eos_evaluators.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

