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
    COM_MillingtonQuirk.cc
    DensityTabular.cc
    EOSDensityEvaluator.cc
    EOSViscosityEvaluator.cc
    H2O_Density.cc
    H2O_DensityFEHM.cc
    H2O_ThermalConductivity.cc
    H2O_SaturatedVaporPressure.cc
    H2O_Viscosity.cc
    H2O_ViscosityFEHM.cc
    Henry_Solubility.cc
    IdealGas_Density.cc
    IsobaricEOSEvaluator.cc
    LookupTable.cc
    MolarFractionGasEvaluator.cc
    VaporInGas_Density.cc
    VaporInGas_Diffusion.cc
    VaporLiquid_Tabular.cc
    ViscosityConstant.cc
    ViscosityTabular.cc
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
    COM_MillingtonQuirk.hh
    DensityTabular.hh
    EOSDensityEvaluator.hh
    EOSViscosityEvaluator.hh
    EOS_Density.hh
    EOS_Diffusion.hh
    EOS_SaturatedVaporPressure.hh
    H2O_ThermalConductivity.hh
    EOS_Utils.hh
    EOS_Viscosity.hh
    eos_reg.hh
    EOSFactory.hh
    H2O_Density.hh
    H2O_DensityFEHM.hh
    H2O_SaturatedVaporPressure.hh
    H2O_Viscosity.hh
    H2O_ViscosityFEHM.hh
    Henry_Solubility.hh
    IdealGas_Density.hh
    IsobaricEOSEvaluator.hh
    LookupTable.hh
    MolarFractionGasEvaluator.hh
    VaporInGas_Density.hh
    VaporInGas_Diffusion.hh
    ViscosityConstant.hh
    ViscosityTabular.hh
   )

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

    set(amanzi_libs eos atk operators)

    # Test: Multiscale
    add_amanzi_test(eos_evaluators eos_evaluators
                    KIND unit
                    SOURCE test/Main.cc test/eos_evaluators.cc
                    LINK_LIBS ${amanzi_libs} ${UnitTest_LIBRARIES})
endif()

