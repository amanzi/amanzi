# -*- mode: cmake -*-

#
#  Amanzi
#   Input examples
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)


set(input_example_src input_example.cpp)
add_executable(example.exe ${input_example_src})
target_link_libraries(example.exe ${Teuchos_LIBRARIES})

add_executable(teuchos_example.exe teuchos_example.cpp)
target_link_libraries(teuchos_example.exe ${Teuchos_LIBRARIES})
