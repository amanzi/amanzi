#set -x
SED=sed
SED_ARGS="-i "

echo "Running sed on $1"

${SED} ${SED_ARGS} -E 's|AmanziMesh::CELL|AmanziMesh::Entity_kind::CELL|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::FACE|AmanziMesh::Entity_kind::FACE|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::NODE|AmanziMesh::Entity_kind::NODE|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::EDGE|AmanziMesh::Entity_kind::EDGE|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::BOUNDARY_FACE|AmanziMesh::Entity_kind::BOUNDARY_FACE|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::PRISM|AmanziMesh::Cell_type::PRISM|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::POLYHED|AmanziMesh::Cell_type::POLYHED|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::HEX|AmanziMesh::Cell_type::HEX|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::PYRAMID|AmanziMesh::Cell_type::PYRAMID|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::QUAD|AmanziMesh::Cell_type::QUAD|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::POLYGON|AmanziMesh::Cell_type::POLYGON|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::TET|AmanziMesh::Cell_type::TET|g' "$1"
${SED} ${SED_ARGS} -E 's|AmanziMesh::TRI|AmanziMesh::Cell_type::TRI|g' "$1"

${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)face_get_cells\(([^,\)]*),([^,\)]*),[ ]*&([^,\)]*)\)|\1\5 = \2getFaceCells\(\3,\4\)|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|(^.*) (\S*)cell_get_faces\(([^,\)]*),[ ]*&([^,\)]*)\)|\1 \4 = \2getCellFaces\(\3\)|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)node_get_cells\(([^,\)]*),([^,\)]*),[ ]*&([^,\)]*)\)|\1\5 = \2getNodeCells\(\3,\4\)|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)edge_get_cells\(([^,\)]*),([^,\)]*),[ ]*&([^,\)]*)\)|\1\5 = \2getEdgeCells\(\3,\4\)|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)cell_get_edges\(([^,\)]*),[ ]*&([^,\)]*)\)|\1\4 = \2getCellEdges\(\3\)|g' "$1" # Need change for type 

${SED} ${SED_ARGS} -E 's|get_set_size\(|getSetSize\(|g' "$1"
${SED} ${SED_ARGS} -E 's|get_set_entities\(|getSetEntities\(|g' "$1"
${SED} ${SED_ARGS} -E 's|get_set_entities_and_vofs\(|getSetEntitiesAndVolumeFractions\(|g' "$1"
${SED} ${SED_ARGS} -E 's|get_comm\(|getComm\(|g' "$1"
${SED} ${SED_ARGS} -E 's|valid_set_name|isValidSetName|g' "$1"
${SED} ${SED_ARGS} -E 's|([^_])cell_volume\(|\1getCellVolume\(|g' "$1"
${SED} ${SED_ARGS} -E 's|face_area\(|getFaceArea\(|g' "$1"
${SED} ${SED_ARGS} -E 's|face_get_cells|getFaceCells|g' "$1"
${SED} ${SED_ARGS} -E 's|\"MeshLight.hh\"|\"Mesh.hh\"|g' "$1"
${SED} ${SED_ARGS} -E 's|MeshLight|Mesh|g' "$1"

${SED} ${SED_ARGS} -E 's|space_dimension\(|getSpaceDimension\(|g' "$1"

# 1 change in one line
${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)cell_get_nodes\(([^,\)]*),[ ]*&([^,\)]*)\)|\1\4 = \2getCellNodes\(\3\)|g' "$1" # Need change for type 
# 2 change multiline 
${SED} ${SED_ARGS} -E 's|cell_get_nodes\(|getCellNodes\(|g' "$1"
${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)node_get_coordinates\(([^,\)]*),[ ]*&([^,\)]*)\)|\1\4 = \2getNodeCoordinate\(\3\)|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|node_get_coordinates\(|getNodeCoordinate\(|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|(^[ \t]*)(.*)face_get_nodes\(([^,\)]*),[ ]*&([^,\)]*)\)|\1\4 = \2getFaceNodes\(\3\)|g' "$1" # Need change for type 
${SED} ${SED_ARGS} -E 's|face_get_nodes\(|getFaceNodes\(|g' "$1" 

${SED} ${SED_ARGS} -E 's|cell_get_faces\(|getCellFaces\(|g' "$1"
${SED} ${SED_ARGS} -E 's|node_get_cells\(|getNodeCells\(|g' "$1"
${SED} ${SED_ARGS} -E 's|edge_get_cells\(|getEdgeCells\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_centroid\(|getCellCentroid\(|g' "$1"
${SED} ${SED_ARGS} -E 's|face_centroid\(|getFaceCentroid\(|g' "$1"
${SED} ${SED_ARGS} -E 's|edge_centroid\(|getEdgeCentroid\(|g' "$1"
${SED} ${SED_ARGS} -E 's|edge_get_nodes\(|getEdgeNodes\(|g' "$1"

${SED} ${SED_ARGS} -E 's|face_normal\(|getFaceNormal\(|g' "$1"
${SED} ${SED_ARGS} -E 's|getFaceNormal\(([^,\)]*),([^,\)]*),([^,\)]*),([^,\)]*)\)|getFaceNormal\(\1,\3,\4\)|g' "$1"
${SED} ${SED_ARGS} -E 's|getFaceNormal\(([^,\)]*),([^,\)]*),([^,\)]*)\)|getFaceNormal\(\1,\3\)|g' "$1"


${SED} ${SED_ARGS} -E 's|cell_get_face_dirs\(|getCellFacesAndDirections\(|g' "$1"
${SED} ${SED_ARGS} -E 's|face_get_edges_and_dirs\(|getFaceEdgesAndDirections\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_edges\(|getCellEdges\(|g' "$1"
${SED} ${SED_ARGS} -E 's|edge_vector\(|getEdgeVector\(|g' "$1"
${SED} ${SED_ARGS} -E 's|edge_length\(|getEdgeLength\(|g' "$1"
${SED} ${SED_ARGS} -E 's|face_to_cell_edge_map\(|getFaceCellEdgeMap\(|g' "$1" 
${SED} ${SED_ARGS} -E 's|cell_get_num_faces\(|getCellNumFaces\(|g' "$1" 
${SED} ${SED_ARGS} -E 's|num_entities\(|getNumEntities\(|g' "$1" 
${SED} ${SED_ARGS} -E 's|edge_get_ho_nodes\(|getEdgeHOCoordinates\(|g' "$1" 
${SED} ${SED_ARGS} -E 's|face_get_ho_nodes\(|getFaceHOCoordinates\(|g' "$1" 
${SED} ${SED_ARGS} -E 's|cell_get_faces_and_dirs\(|getCellFacesAndDirections\(|g' "$1"

${SED} ${SED_ARGS} -E 's|->importer\(|->getImporter\(|g' "$1"
${SED} ${SED_ARGS} -E 's|\.importer\(|\.getImporter\(|g' "$1"
${SED} ${SED_ARGS} -E 's|exterior_face_importer\(|getBoundaryFaceImporter\(|g' "$1"

${SED} ${SED_ARGS} -E 's|([^_])face_map\(([^,\)]*)\)|\1getMap\(AmanziMesh::Entity_kind::FACE,\2\)|g' "$1"
${SED} ${SED_ARGS} -E 's|([^_])node_map\(([^,\)]*)\)|\1getMap\(AmanziMesh::Entity_kind::NODE,\2\)|g' "$1"
${SED} ${SED_ARGS} -E 's|([^_])edge_map\(([^,\)]*)\)|\1getMap\(AmanziMesh::Entity_kind::EDGE,\2\)|g' "$1"
${SED} ${SED_ARGS} -E 's|([^_])cell_map\(([^,\)]*)\)|\1getMap\(AmanziMesh::Entity_kind::CELL,\2\)|g' "$1"
${SED} ${SED_ARGS} -E 's|exterior_face_map\(([^,\)]*)\)|getMap\(AmanziMesh::Entity_kind::BOUNDARY_FACE,\1\)|g' "$1"
${SED} ${SED_ARGS} -E 's|->map\(|->getMap\(|g' "$1"
${SED} ${SED_ARGS} -E 's|\.map\(|\.getMap\(|g' "$1"

${SED} ${SED_ARGS} -E 's|entity_kind_string\(|to_string\(|g' "$1"
${SED} ${SED_ARGS} -E 's|=.*cell_type_to_name\(|= to_string\(|g' "$1"
${SED} ${SED_ARGS} -E 's|manifold_dimension\(|getManifoldDimension\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_max_nodes\(|getCellMaxNodes\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_max_faces\(|getCellMaxFaces\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_max_edges\(|getCellMaxEdges\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_type\(|getCellType\(|g' "$1"
${SED} ${SED_ARGS} -E 's|is_logical\(|isLogical\(|g' "$1"
${SED} ${SED_ARGS} -E 's|vis_mesh\(|getVisMesh\(|g' "$1"

${SED} ${SED_ARGS} -E 's|get_referencing_parent\(|getReferencingParent\(|g' "$1"
${SED} ${SED_ARGS} -E 's|entity_kind\(|createEntityKind\(|g' "$1"
${SED} ${SED_ARGS} -E 's|DoImport\(|doImport\(|g' "$1"
${SED} ${SED_ARGS} -E 's|geometric_model\(|getGeometricModel\(|g' "$1"
${SED} ${SED_ARGS} -E 's|build_columns\(|buildColumns\(|g' "$1"
${SED} ${SED_ARGS} -E 's|num_columns\(|columns\(|g' "$1"
${SED} ${SED_ARGS} -E 's|parameter_list\(|getParameterList\(|g' "$1"
${SED} ${SED_ARGS} -E 's|valid_edges\(|hasEdges\(|g' "$1"
${SED} ${SED_ARGS} -E 's|entity_get_parent\(|getEntityParent\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_faces_and_bisectors\(|getCellFacesAndBisectors\(|g' "$1"
${SED} ${SED_ARGS} -E 's|DoExport\(|doExport\(|g' "$1"
${SED} ${SED_ARGS} -E 's|parent\(|getParentMesh\(|g' "$1"
${SED} ${SED_ARGS} -E 's|cell_get_face_adj_cells\(|AmanziMesh::MeshAlgorithms::getCellFaceAdjacentCells\(|g' "$1"
${SED} ${SED_ARGS} -E 's|get_indexing_parent\(|getIndexingParent\(|g' "$1"


git -P diff --shortstat $1