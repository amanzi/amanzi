# -*- mode: cmake -*-

#
#  Amanzi
#    Mesh 
#
#

# Amanzi module, include files found in AMANZI_MODULE_PATH
include(PrintVariable)

#
# Define a project name
# After this command the following varaibles are defined
#   MESH_SOURCE_DIR
#   MESH_BINARY_DIR
# Other projects (subdirectories) can reference this directory
# through these variables.
project(MESH)

# Amanzi include directories
include_directories(${DBC_SOURCE_DIR})
include_directories(${GEOMETRY_SOURCE_DIR})

# External (TPL) directories
include_directories(${Epetra_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})

#
# Library: mesh
#
add_library(mesh Mesh.cc GenerationSpec.cc)
add_install_library(mesh)
target_link_libraries(mesh geometry)



#
# Library: cell_geometry (WILL GO AWAY WHEN EVERYONE SWITCHES TO THE NEW API)
#
add_library(cell_geometry cell_geometry.cc)
add_install_library(cell_geometry)

#
# Library: mesh_audit
#
add_library(mesh_audit MeshAudit.cc MeshAuditOld.cc)
add_install_library(mesh_audit)
target_link_libraries(mesh_audit mesh error_handling cell_geometry ${Epetra_LIBRARIES})

#
# Install Header files
#
file(GLOB mesh_inc_files "*.hh")
add_install_include_file(${mesh_inc_files})



# Derived and support classes

add_subdirectory(mesh_data)
add_subdirectory(mesh_simple)
add_subdirectory(exodus_reader)

# Mesh Frameworks

# STK (Trilinos Package)
if (ENABLE_STK_Mesh)
    add_subdirectory(mesh_stk)
endif(ENABLE_STK_Mesh)    

# MOAB Mesh framework
if (ENABLE_MOAB_Mesh)
    add_subdirectory(mesh_moab)
endif(ENABLE_MOAB_Mesh)    

# MSTK Mesh framework
if (ENABLE_MSTK_Mesh)
    add_subdirectory(mesh_mstk)
endif(ENABLE_MSTK_Mesh)    


# Factory class

add_subdirectory(mesh_factory)

