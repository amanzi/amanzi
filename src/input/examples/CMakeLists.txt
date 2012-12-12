# -*- mode: cmake -*-

#
#  Amanzi
#   Input examples
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

# Add the input project directory to the include list
include_directories(${INPUT_SOURCE_DIR})

set(input_example_src input_example.cc)
add_executable(input_example.exe ${input_example_src})
#target_link_libraries(input_example.exe input input_test ${Teuchos_LIBRARIES})
target_link_libraries(input_example.exe input ${Teuchos_LIBRARIES})

add_executable(teuchos_example.exe teuchos_example.cpp)
target_link_libraries(teuchos_example.exe ${Teuchos_LIBRARIES})
