function go () {
    echo "Running go on $1"
    # simple renames
    gsed -i 's|\.num_entities(|.getNumEntities(|g' "$1"
    gsed -i 's|->num_entities(|->getNumEntities(|g' "$1"
    gsed -i 's|entity_get_ptype(|getEntityPtype(|g' "$1"
    gsed -i 's|entity_get_parent(|getEntityParent(|g' "$1"
    gsed -i 's|node_set_coordinates|setNodeCoordinate|g' "$1"

    # mesh.node_get_coordinate(n, &node) --> node = mesh.getNodeCoordinate(n)
    # for node, face, cell
    gsed -iE 's|\([^ ]*\)node_get_coordinates(\([^,]*\),[ ]*&\(.*\));|\3 = \1getNodeCoordinate(\2);|g' "$1"
    gsed -iE 's|\([^ ]*\)cell_get_coordinates(\([^,]*\),[ ]*&\(.*\));|\3 = \1getCellCoordinates(\2);|g' "$1"
    gsed -iE 's|\([^ ]*\)face_get_coordinates(\([^,]*\),[ ]*&\(.*\));|\3 = \1getFaceCoordinates(\2);|g' "$1"
    
    # deal with both removing recompute option and case without that option
    gsed -i 's@\(\.\|->\)cell_centroid(\([^,]*\)[ ,]*\(true\|false\|\));@\1getCellCentroid(\2)@g' "$1"
    gsed -i 's@\(\.\|->\)cell_volume(\([^,]*\)[ ,]*\(true\|false\|\));@\1getCellVolume(\2)@g' "$1"
    gsed -i 's@\(\.\|->\)face_centroid(\([^,]*\)[ ,]*\(true\|false\|\));@\1getFaceCentroid(\2)@g' "$1"
    gsed -i 's@\(\.\|->\)face_area(\([^,]*\)[ ,]*\(true\|false\|\));@\1getFaceArea(\2)@g' "$1"
    gsed -i 's@\(\.\|->\)face_normal(\([^,]*\)[ ,]*\(true\|false\|\));@\1getFaceNormal(\2)@g' "$1"
    gsed -i 's@\(\.\|->\)edge_vector(\([^,]*\)[ ,]*\(true\|false\|\));@\1getEdgeVector(\2)@g' "$1"
    gsed -i 's@\(\.\|->\)edge_length(\([^,]*\)[ ,]*\(true\|false\|\));@\1getEdgeLength(\2)@g' "$1"

    gsed -i 's@\.edge_centroid(@.getEdgeCentroid(@g' "$1"
    gsed -i 's@->edge_centroid(@.getEdgeCentroid(@g' "$1"
    
    # remove ordered option, remove reference
    # cell_get_faces_and_dirs(c, &faces, &dirs, true)
    #   --> getCellFacesAndDirs(c, faces, &dirs)
    gsed -iE 's@cell_get_faces_and_dirs(\([^&]*\)&\(.*\),[ ]*\(true\|false\)[ ]*)@getCellFacesAndDirs(\1\2)@g' "$1"
    gsed -iE 's@cell_get_faces_and_dirs(\([^&]*\)&\(.*\))@getCellFacesAndDirs(\1\2)@g' "$1"

    gsed -iE 's@cell_get_faces_and_bisectors(\([^&]*\)&\(.*\),[ ]*\(true\|false\)[ ]*)@getCellFacesAndBisectors(\1\2)@g' "$1"
    gsed -iE 's@cell_get_faces_and_bisectors(\([^&]*\)&\(.*\))@getCellFacesAndBisectors(\1\2)@g' "$1"
    
    gsed -iE 's@face_get_edges_and_dirs(\([^&]*\)&\(.*\),[ ]*\(true\|false\)[ ]*)@getFaceEdgesAndDirs(\1\2)@g' "$1"
    gsed -iE 's@face_get_edges_and_dirs(\([^&]*\)&\(.*\))@getFaceEdgesAndDirs(\1\2)@g' "$1"

    # old style cell_get_faces
    gsed -iE 's|cell_get_faces(\([^&]*\)&\(.*\))|getCellFaces(\1\2)|g' "$1"
    # new style cell get faces
    gsed -iE 's|=\(.*\)cell_get_faces|=\1getCellFaces|g' "$1"
    
    gsed -iE 's|cell_get_edges(\([^&]*\)&\(.*\))|getCellEdges(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)cell_get_edges|=\1getCellEdges|g' "$1"
    gsed -iE 's|cell_get_nodes(\([^&]*\)&\(.*\))|getCellNodes(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)cell_get_nodes|=\1getCellNodes|g' "$1"
    gsed -iE 's|face_get_edges(\([^&]*\)&\(.*\))|getFaceEdges(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)face_get_edges|=\1getFaceEdges|g' "$1"
    gsed -iE 's|face_get_nodes(\([^&]*\)&\(.*\))|getFaceNodes(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)face_get_nodes|=\1getFaceNodes|g' "$1"
    gsed -iE 's|edge_get_nodes(\([^&]*\)&\(.*\))|getEdgeNodes(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)edge_get_nodes|=\1getEdgeNodes|g' "$1"

    gsed -iE 's|face_get_cells(\([^&]*\)&\(.*\))|getFaceCells(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)face_get_cells|=\1getFaceCells|g' "$1"
    gsed -iE 's|edge_get_cells(\([^&]*\)&\(.*\))|getEdgeCells(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)edge_get_cells|=\1getEdgeCells|g' "$1"
    gsed -iE 's|node_get_cells(\([^&]*\)&\(.*\))|getNodeCells(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)node_get_cells|=\1getNodeCells|g' "$1"
    gsed -iE 's|edge_get_faces(\([^&]*\)&\(.*\))|getEdgeFaces(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)edge_get_faces|=\1getEdgeFaces|g' "$1"
    gsed -iE 's|node_get_faces(\([^&]*\)&\(.*\))|getNodeFaces(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)node_get_faces|=\1getNodeFaces|g' "$1"
    gsed -iE 's|node_get_edges(\([^&]*\)&\(.*\))|getNodeEdges(\1\2)|g' "$1"
    gsed -iE 's|=\(.*\)node_get_edges|=\1getNodeEdges|g' "$1"

    gsed -i "s/AmanziMesh::CELL/AmanziMesh::Entity_kind::CELL/g" "$1"
    gsed -i "s/AmanziMesh::FACE/AmanziMesh::Entity_kind::FACE/g" "$1"
    gsed -i "s/AmanziMesh::EDGE/AmanziMesh::Entity_kind::EDGE/g" "$1"
    gsed -i "s/AmanziMesh::NODE/AmanziMesh::Entity_kind::NODE/g" "$1"
    gsed -i "s/AmanziMesh::BOUNDARY_FACE/AmanziMesh::Entity_kind::BOUNDARY_FACE/g" "$1"
    gsed -i "s/Mesh::TRI/Mesh::Cell_type::TRI/g" "$1"
    gsed -i "s/Mesh::QUAD/Mesh::Cell_type::QUAD/g" "$1"
    gsed -i "s/Mesh::POLYGON/Mesh::Cell_type::POLYGON/g" "$1"
    gsed -i "s/Mesh::TET/Mesh::Cell_type::TET/g" "$1"
    gsed -i "s/Mesh::PRISM/Mesh::Cell_type::PRISM/g" "$1"
    gsed -i "s/Mesh::PYRAMID/Mesh::Cell_type::PYRAMID/g" "$1"
    gsed -i "s/Mesh::HEX/Mesh::Cell_type::HEX/g" "$1"
    gsed -i "s/Mesh::POLYHED/Mesh::Cell_type::POLYHED/g" "$1"

    gsed -i "s^manifold_dimension^get_manifold_dimension^g" "$1"
    gsed -i "s^space_dimension^get_space_dimension^g" "$1"

    gsed -i 's^"MeshFramework\.hh"^"MeshFrameworkTraits.hh"^g' "$1"
    gsed -i 's^"Mesh\.hh"^"MeshFramework.hh"^g' "$1"

    # guess some mesh GID calls -- cannot blanket replace because Map.GID()
    gsed -i 's^mesh->GID^mesh->getEntityGID^g' "$1"
    gsed -i 's^mesh\.GID^mesh.getEntityGID^g' "$1"
    gsed -i 's^mesh_->GID^mesh_->getEntityGID^g' "$1"
    gsed -i 's^mesh_\.GID^mesh_.getEntityGID^g' "$1"
}
