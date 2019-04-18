# -*- mode: cmake -*-
#
#  ATS
#    Constitutive relations for flow
#
include_directories(${Amanzi_TPL_MSTK_INCLUDE_DIRS})
include_directories(${ATS_SOURCE_DIR}/src/factory)
add_definitions("-DMSTK_HAVE_MPI")

list(APPEND subdirs wrm porosity overland_conductivity elevation water_content sources thaw_depth)

set(sources "")

# include all subdirectories
foreach(lcv IN LISTS subdirs)
  include_directories($lcv)
endforeach(lcv)

# collect all registrations and register
file(GLOB_RECURSE registrations "./*_reg.hh" )
foreach(reg_lcv IN LISTS registrations)
  register_abs_evaluator_with_factory(HEADERFILE ${reg_lcv} LISTNAME FLOW_RELATIONS)
endforeach(reg_lcv)

generate_evaluators_registration_header(
  HEADERFILE flow_relations_registration.hh
  LISTNAME   FLOW_RELATIONS
  INSTALL    True
  )


# collect all sources
foreach(lcv IN LISTS subdirs)
  file(GLOB subdir_sources "./${lcv}/*.cc")
  set(sources ${sources} ${subdir_sources})
endforeach(lcv)

# make the library
add_library(flow_relations ${sources})

	
if (BUILD_TESTS)
    # Add UnitTest includes
    include_directories(${Amanzi_TPL_UnitTest_INCLUDE_DIRS})
    include_directories("./wrm")

    add_amanzi_test(wrm_vanGenuchten wrm_vanGenuchten
                    KIND unit
                    SOURCE wrm/models/test/main.cc
                           wrm/models/test/test_vanGenuchten.cc
                    LINK_LIBS flow_relations amanzi_error_handling amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES} ${Amanzi_TPL_Trilinos_LIBRARIES})

    add_amanzi_test(wrm_plantChristoffersen wrm_plantChristoffersen
                    KIND unit
                    SOURCE wrm/models/test/main.cc
                           wrm/models/test/test_vanGenuchten.cc
                    LINK_LIBS flow_relations amanzi_error_handling amanzi_state ${Amanzi_TPL_UnitTest_LIBRARIES} ${Amanzi_TPL_Trilinos_LIBRARIES})


endif()
