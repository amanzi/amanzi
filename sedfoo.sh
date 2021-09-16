
function go () {
    SED=gsed
    SED_ARGS=-i''

    echo "Running go on $1"
    # simple renames
    ${SED} ${SED_ARGS} 's|\.num_entities(|.getNumEntities(|g' "$1"
    ${SED} ${SED_ARGS} 's|->num_entities(|->getNumEntities(|g' "$1"
    ${SED} ${SED_ARGS} 's|entity_get_ptype(|getEntityPtype(|g' "$1"
    ${SED} ${SED_ARGS} 's|entity_get_parent(|getEntityParent(|g' "$1"
    ${SED} ${SED_ARGS} 's|node_set_coordinates|setNodeCoordinate|g' "$1"

    ${SED} ${SED_ARGS} 's|\.cell_map(|.getMap(AmanziMesh::Entity_kind::CELL, |g' "$1"
    ${SED} ${SED_ARGS} 's|->cell_map(|->getMap(AmanziMesh::Entity_kind::CELL, |g' "$1"
    ${SED} ${SED_ARGS} 's|\.face_map(|.getMap(AmanziMesh::Entity_kind::FACE, |g' "$1"
    ${SED} ${SED_ARGS} 's|->face_map(|->getMap(AmanziMesh::Entity_kind::FACE, |g' "$1"
    ${SED} ${SED_ARGS} 's|\.edge_map(|.getMap(AmanziMesh::Entity_kind::EDGE, |g' "$1"
    ${SED} ${SED_ARGS} 's|->edge_map(|->getMap(AmanziMesh::Entity_kind::EDGE, |g' "$1"
    ${SED} ${SED_ARGS} 's|\.node_map(|.getMap(AmanziMesh::Entity_kind::NODE, |g' "$1"
    ${SED} ${SED_ARGS} 's|->node_map(|->getMap(AmanziMesh::Entity_kind::NODE, |g' "$1"
    ${SED} ${SED_ARGS} 's|\.exterior_face_map(|.getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, |g' "$1"
    ${SED} ${SED_ARGS} 's|->exterior_face_map(|->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE, |g' "$1"
    ${SED} ${SED_ARGS} 's|\.exterior_node_map(|.getMap(AmanziMesh::Entity_kind::BOUNDARY_NODE, |g' "$1"
    ${SED} ${SED_ARGS} 's|->exterior_node_map(|->getMap(AmanziMesh::Entity_kind::BOUNDARY_NODE, |g' "$1"
    ${SED} ${SED_ARGS} 's|->get_comm(|->getComm(|g' "$1"
    ${SED} ${SED_ARGS} 's|\.get_comm(|.getComm(|g' "$1"
    ${SED} ${SED_ARGS} 's|->map(|->getMap(|g' "$1"
    ${SED} ${SED_ARGS} 's|\.map(|.getMap(|g' "$1"

    # mesh.node_get_coordinate(n, &node) --> node = mesh.getNodeCoordinate(n)
    ${SED} ${SED_ARGS} -E 's@([^ ]*)node_get_coordinate\(([^,]),[ ]*&(.*)\)@\3 = \1getNodeCoordinate(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)cell_get_coordinates\(([^,]),[ ]*&(.*)\)@\3 = \1getCellCoordinates(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)face_get_coordinates\(([^,]),[ ]*&(.*)\)@\3 = \1getFaceCoordinates(\2)@g' "$1"

    # deal with both removing recompute option and case without that option
    #
    # mesh->cell_centroid(c, true) --> mesh->getCellCentroid(c)
    # mesh->cell_centroid(c) --> mesh->getCellCentroid(c)

    ${SED} ${SED_ARGS} -E 's@([^ ]*)cell_centroid\(([^,]*)[, ]*(true|false|)\)@\1getCellCentroid(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)cell_volume\(([^,]*)[, ]*(true|false|)\)@\1getCellVolume(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)face_centroid\(([^,]*)[, ]*(true|false|)\)@\1getFaceCentroid(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)face_area\(([^,]*)[, ]*(true|false|)\)@\1getFaceArea(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)edge_length\(([^,]*)[, ]*(true|false|)\)@\1getEdgeLength(\2)@g' "$1"

    # above, but also deals with extra arguments afterwards -- first
    # an id, then an orientation.  We deal with the case of with the ID first
    ${SED} ${SED_ARGS} -E 's@([^ ]*)edge_vector\(([^,]*)[, ]*(true|false),(.*)\)@\1getEdgeVector(\2, \4)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)face_normal\(([^,]*)[, ]*(true|false),(.*)\)@\1getFaceNormal(\2, \4)@g' "$1"
    # then the case without ID or orientation
    ${SED} ${SED_ARGS} -E 's@([^ ]*)edge_vector\(([^,]*)[, ]*(true|false|)\)@\1getEdgeVector(\2)@g' "$1"
    ${SED} ${SED_ARGS} -E 's@([^ ]*)face_normal\(([^,]*)[, ]*(true|false|)\)@\1getFaceNormal(\2)@g' "$1"

    
    
    ${SED} ${SED_ARGS} 's@\.edge_centroid(@.getEdgeCentroid(@g' "$1"
    ${SED} ${SED_ARGS} 's@->edge_centroid(@.getEdgeCentroid(@g' "$1"
    
    # remove ordered option, remove reference
    # cell_get_faces_and_dirs(c, &faces, &dirs, true)
    #   --> getCellFacesAndDirs(c, faces, &dirs)
    ${SED} ${SED_ARGS} 's@cell_get_faces_and_dirs(\([^&]*\)&\(.*\),[ ]*\(true\|false\)[ ]*)@getCellFacesAndDirs(\1\2)@g' "$1"
    ${SED} ${SED_ARGS} 's@cell_get_faces_and_dirs(\([^&]*\)&\(.*\))@getCellFacesAndDirs(\1\2)@g' "$1"

    ${SED} ${SED_ARGS} 's@cell_get_faces_and_bisectors(\([^&]*\)&\(.*\),[ ]*\(true\|false\)[ ]*)@getCellFacesAndBisectors(\1\2)@g' "$1"
    ${SED} ${SED_ARGS} 's@cell_get_faces_and_bisectors(\([^&]*\)&\(.*\))@getCellFacesAndBisectors(\1\2)@g' "$1"
    
    ${SED} ${SED_ARGS} 's@face_get_edges_and_dirs(\([^&]*\)&\(.*\),[ ]*\(true\|false\)[ ]*)@getFaceEdgesAndDirs(\1\2)@g' "$1"
    ${SED} ${SED_ARGS} 's@face_get_edges_and_dirs(\([^&]*\)&\(.*\))@getFaceEdgesAndDirs(\1\2)@g' "$1"

    # old style cell_get_faces
    ${SED} ${SED_ARGS} 's|cell_get_faces(\([^&]*\)&\(.*\))|getCellFaces(\1\2)|g' "$1"
    # new style cell get faces
    ${SED} ${SED_ARGS} 's|=\(.*\)cell_get_faces|=\1getCellFaces|g' "$1"
    
    ${SED} ${SED_ARGS} 's|cell_get_edges(\([^&]*\)&\(.*\))|getCellEdges(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)cell_get_edges|=\1getCellEdges|g' "$1"
    ${SED} ${SED_ARGS} 's|cell_get_nodes(\([^&]*\)&\(.*\))|getCellNodes(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)cell_get_nodes|=\1getCellNodes|g' "$1"
    ${SED} ${SED_ARGS} 's|face_get_edges(\([^&]*\)&\(.*\))|getFaceEdges(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)face_get_edges|=\1getFaceEdges|g' "$1"
    ${SED} ${SED_ARGS} 's|face_get_nodes(\([^&]*\)&\(.*\))|getFaceNodes(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)face_get_nodes|=\1getFaceNodes|g' "$1"
    ${SED} ${SED_ARGS} 's|edge_get_nodes(\([^&]*\)&\(.*\))|getEdgeNodes(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)edge_get_nodes|=\1getEdgeNodes|g' "$1"

    ${SED} ${SED_ARGS} 's|face_get_cells(\([^&]*\)&\(.*\))|getFaceCells(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)face_get_cells|=\1getFaceCells|g' "$1"
    ${SED} ${SED_ARGS} 's|edge_get_cells(\([^&]*\)&\(.*\))|getEdgeCells(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)edge_get_cells|=\1getEdgeCells|g' "$1"
    ${SED} ${SED_ARGS} 's|node_get_cells(\([^&]*\)&\(.*\))|getNodeCells(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)node_get_cells|=\1getNodeCells|g' "$1"
    ${SED} ${SED_ARGS} 's|edge_get_faces(\([^&]*\)&\(.*\))|getEdgeFaces(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)edge_get_faces|=\1getEdgeFaces|g' "$1"
    ${SED} ${SED_ARGS} 's|node_get_faces(\([^&]*\)&\(.*\))|getNodeFaces(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)node_get_faces|=\1getNodeFaces|g' "$1"
    ${SED} ${SED_ARGS} 's|node_get_edges(\([^&]*\)&\(.*\))|getNodeEdges(\1\2)|g' "$1"
    ${SED} ${SED_ARGS} 's|=\(.*\)node_get_edges|=\1getNodeEdges|g' "$1"

    ${SED} ${SED_ARGS} "s/AmanziMesh::CELL/AmanziMesh::Entity_kind::CELL/g" "$1"
    ${SED} ${SED_ARGS} "s/AmanziMesh::FACE/AmanziMesh::Entity_kind::FACE/g" "$1"
    ${SED} ${SED_ARGS} "s/AmanziMesh::EDGE/AmanziMesh::Entity_kind::EDGE/g" "$1"
    ${SED} ${SED_ARGS} "s/AmanziMesh::NODE/AmanziMesh::Entity_kind::NODE/g" "$1"
    ${SED} ${SED_ARGS} "s/AmanziMesh::BOUNDARY_FACE/AmanziMesh::Entity_kind::BOUNDARY_FACE/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::TRI/Mesh::Cell_type::TRI/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::QUAD/Mesh::Cell_type::QUAD/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::POLYGON/Mesh::Cell_type::POLYGON/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::TET/Mesh::Cell_type::TET/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::PRISM/Mesh::Cell_type::PRISM/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::PYRAMID/Mesh::Cell_type::PYRAMID/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::HEX/Mesh::Cell_type::HEX/g" "$1"
    ${SED} ${SED_ARGS} "s/Mesh::POLYHED/Mesh::Cell_type::POLYHED/g" "$1"

    ${SED} ${SED_ARGS} "s^manifold_dimension^getManifoldDimension^g" "$1"
    ${SED} ${SED_ARGS} "s^space_dimension^getSpaceDimension^g" "$1"

    ${SED} ${SED_ARGS} 's^"MeshFramework\.hh"^"MeshFrameworkTraits.hh"^g' "$1"
    ${SED} ${SED_ARGS} 's^"Mesh\.hh"^"MeshFramework.hh"^g' "$1"

    # guess some mesh GID calls -- cannot blanket replace because Map.GID()
    ${SED} ${SED_ARGS} 's^mesh->GID^mesh->getEntityGID^g' "$1"
    ${SED} ${SED_ARGS} 's^mesh\.GID^mesh.getEntityGID^g' "$1"
    ${SED} ${SED_ARGS} 's^mesh_->GID^mesh_->getEntityGID^g' "$1"
    ${SED} ${SED_ARGS} 's^mesh_\.GID^mesh_.getEntityGID^g' "$1"

}
