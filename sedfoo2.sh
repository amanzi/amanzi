# this fixed some intermediate changes of mesh method names
function go2 () {
    SED=gsed
    SED_ARGS=-i''

    echo "Running go2 on $1"

    ${SED} ${SED_ARGS} 's|et_comm|etComm|g' "$1"
    ${SED} ${SED_ARGS} 's|et_parameter_list|etParameterList|g' "$1"
    ${SED} ${SED_ARGS} 's|et_parameter_list|etParameterList|g' "$1"
    ${SED} ${SED_ARGS} 's|et_verbose_object|etVerboseObject|g' "$1"
    ${SED} ${SED_ARGS} 's|et_geometric_model|etGeometricModel|g' "$1"
    ${SED} ${SED_ARGS} 's|et_parent|etParentMesh|g' "$1"
    ${SED} ${SED_ARGS} 's|et_vis_mesh|etVisMesh|g' "$1"
    ${SED} ${SED_ARGS} 's|has_edges|hasEdges|g' "$1"
    ${SED} ${SED_ARGS} 's|is_ordered|isOrdered|g' "$1"
    ${SED} ${SED_ARGS} 's|is_deformable|isDeformable|g' "$1"
    ${SED} ${SED_ARGS} 's|get_importer|getImporter|g' "$1"

    ${SED} ${SED_ARGS} 's^"MeshFramework\.hh"^"MeshFrameworkTraits.hh"^g' "$1"
    ${SED} ${SED_ARGS} 's^"Mesh\.hh"^"MeshFramework.hh"^g' "$1"

}
    
