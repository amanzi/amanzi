#
# ATS
#
# NOTE: ATS is not standalone code, and this is not a stand-alone
#       CMakeLists.txt.  Instead, it must be built as a subrepo of
#       Amanzi.  See ATS/INSTALL.md

message(STATUS "\n-------------------------------------------")
message(STATUS "\n-- CMake: Configuring ATS build/install.\n--")
message(STATUS "----------------------------------------")
add_subdirectory(src)

