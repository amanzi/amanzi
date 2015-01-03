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
set(eos_src_files eos_factory.cc 
    eos_evaluator.cc
   )

add_amanzi_library(eos
                   SOURCE ${eos_src_files} HEADERS ${eos_inc_files}
                   LINK_LIBS atk)

#
# Install Targets
#
set(eos_inc_files 
    eos.hh
    eos_factory.hh
    eos_evaluator.hh)

add_install_include_file(${eos_inc_files})

if (BUILD_TESTS) 
endif()

