# -*- mode: cmake -*-
include(TestManager)
include(LibraryManager)

include(PreprocessBoxLibFortran)

# Set locations for amr-related source files, make available outside this dir afterward
set(PMAMR_REL_INCLUDE_DIRS mpc_pk problem_setup simulation_driver utils materials mpc_pk/transport/diffusion)
set(PMAMR_INCLUDE_DIRS)
foreach(d ${PMAMR_REL_INCLUDE_DIRS})
  list(APPEND PMAMR_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/${d})
endforeach(d)

foreach(d ${PMAMR_REL_INCLUDE_DIRS})
  add_subdirectory(${d})
endforeach(d)

add_library(structured_amanzi $<TARGET_OBJECTS:boxprob1> $<TARGET_OBJECTS:structutils1> $<TARGET_OBJECTS:structmat1> $<TARGET_OBJECTS:pmamr1> $<TARGET_OBJECTS:boxdriv1> $<TARGET_OBJECTS:strans1>) 
SET_TARGET_PROPERTIES(structured_amanzi PROPERTIES LINKER_LANGUAGE C)
