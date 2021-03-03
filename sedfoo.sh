function msed () {
    gsed -i "s|$1|$2|g" "$3"
}

function go () {
    gsed -i 's|\.num_entities(|.getNumEntities(|g' "$1"
    gsed -i 's|->num_entities(|->getNumEntities(|g' "$1"

    gsed -i 's|entity_get_ptype(|getEntityPtype(|g' "$1"
    #msed GID getEntityGID "$1"
    gsed -i 's|entity_get_parent(|getEntityParent(|g' "$1"

    gsed -iE 's|\([^ ]*\)node_get_coordinates(\([^,]*\),[ ]*&\(.*\));|\3 = \1getNodeCoordinate(\2)|g'
    gsed -i 's|node_set_coordinates|setNodeCoordinate|g' "$1"
    
    # cell_centroid(c, recompute=true) --> getCellCentroid(c)
    gsed -i 's|\.cell_centroid(\(.*\),[ ]*true);|.getCellCentroid(\1)|g' "$1"
    gsed -i 's|->cell_centroid(\(.*\),[ ]*true);|->getCellCentroid(\1)|g' "$1"
    gsed -i 's|\.cell_centroid(\(.*\),[ ]*false);|.getCellCentroid(\1)|g' "$1"
    gsed -i 's|->cell_centroid(\(.*\),[ ]*false);|->getCellCentroid(\1)|g' "$1"
    gsed -i 's|\.cell_centroid(|.getCellCentroid(|g' "$1"
    gsed -i 's|->cell_centroid(|.getCellCentroid(|g' "$1"

    gsed -i 's|\.cell_volume(\(.*\),[ ]*true);|.getCellVolume(\1)|g' "$1"
    gsed -i 's|->cell_volume(\(.*\),[ ]*true);|->getCellVolume(\1)|g' "$1"
    gsed -i 's|\.cell_volume(\(.*\),[ ]*false);|.getCellVolume(\1)|g' "$1"
    gsed -i 's|->cell_volume(\(.*\),[ ]*false);|->getCellVolume(\1)|g' "$1"
    gsed -i 's|\.cell_volume(|.getCellVolume(|g' "$1"
    gsed -i 's|->cell_volume(|.getCellVolume(|g' "$1"
    
    gsed -i 's|\.face_centroid(\(.*\),[ ]*true);|.getFaceCentroid(\1)|g' "$1"
    gsed -i 's|->face_centroid(\(.*\),[ ]*true);|->getFaceCentroid(\1)|g' "$1"
    gsed -i 's|\.face_centroid(\(.*\),[ ]*false);|.getFaceCentroid(\1)|g' "$1"
    gsed -i 's|->face_centroid(\(.*\),[ ]*false);|->getFaceCentroid(\1)|g' "$1"
    gsed -i 's|\.face_centroid(|.getFaceCentroid(|g' "$1"
    gsed -i 's|->face_centroid(|.getFaceCentroid(|g' "$1"

    gsed -i 's|\.face_normal(\(.*\),[ ]*true\(.*\));|.getFaceNormal(\1\2)|g' "$1"
    gsed -i 's|\->face_normal(\(.*\),[ ]*true\(.*\));|->getFaceNormal(\1\2)|g' "$1"
    gsed -i 's|\.face_normal(\(.*\),[ ]*false\(.*\));|.getFaceNormal(\1\2)|g' "$1"
    gsed -i 's|\->face_normal(\(.*\),[ ]*false\(.*\));|->getFaceNormal(\1\2)|g' "$1"
    gsed -i 's|\.face_normal(|.getFaceNormal(|g' "$1"
    gsed -i 's|->face_normal(|.getFaceNormal(|g' "$1"
    

    
    gsed -i 's|\.edge_centroid(|.getEdgeCentroid(|g' "$1"
    gsed -i 's|->edge_centroid(|.getEdgeCentroid(|g' "$1"
    
    gsed -i 's|\.cell_volume(|.getCellVolume(|g' "$1"
    gsed -i 's|->cell_volume(|.getCellVolume(|g' "$1"

    gsed -i 's|\.face_area(|.getFaceArea(|g' "$1"
    gsed -i 's|->face_area(|.getFaceArea(|g' "$1"
    
    gsed -i 's|\.edge_length(|.getEdgeLength(|g' "$1"
    gsed -i 's|->edge_length(|.getEdgeLength(|g' "$1"
    gsed -i 's|\.edge_vector(|.getEdgeVector(|g' "$1"
    gsed -i 's|->edge_vector(|.getEdgeVector(|g' "$1"
    

    gsed -iE 's|cell_get_faces_and_dirs(\([^&]*\)&\(.*\))|getCellFacesAndDirs(\1\2)|g' "$1"
    gsed -iE 's|cell_get_faces_and_bisectors(\([^&]*\)&\(.*\))|getCellFacesAndBisectors(\1\2)|g' "$1"
    gsed -iE 's|cell_get_faces(\([^&]*\)&\(.*\))|getCellFaces(\1\2)|g' "$1"


    msed cell_get_faces getCellFaces "$1"
    msed cell_get_edges getCellEdges "$1"
    msed cell_get_nodes getCellNodes "$1"
    msed cell_get_2D_edges_and_dirs getCell2DEdgesAndDirs "$1"
    msed face_get_edges_and_dirs getFaceEdgesAndDirs "$1"
    msed face_get_edges getFaceEdges "$1"
    msed face_get_edges getFaceEdges "$1"
    msed face_get_nodes getFaceNodes "$1"
    msed edge_get_nodes getEdgeNodes "$1"
    msed face_get_cells getFaceCells "$1"
    msed edge_get_cells getEdgeCells "$1"
    msed edge_get_faces getEdgeFaces "$1"
    msed node_get_cells getNodeCells "$1"
    msed node_get_faces getNodeFaces "$1"
    msed node_get_edges getNodeEdges "$1"
    msed get_set_entities_and_vofs getSetEntities "$1"
    msed get_set_entities getSetEntities "$1"

    gsed -i "s/AmanziMesh::CELL/AmanziMesh::Entity_kind::CELL/g" "$1"
    gsed -i "s/AmanziMesh::FACE/AmanziMesh::Entity_kind::FACE/g" "$1"
    gsed -i "s/AmanziMesh::EDGE/AmanziMesh::Entity_kind::EDGE/g" "$1"
    gsed -i "s/AmanziMesh::NODE/AmanziMesh::Entity_kind::NODE/g" "$1"
    gsed -i "s/AmanziMesh::BOUNDARY_FACE/AmanziMesh::Entity_kind::BOUNDARY_FACE/g" "$1"

    gsed -i "s^Mesh\.hh^MeshFramework.hh^g" "$1"
    gsed -i "s^MeshFactory\.hh^MeshFrameworkFactory.hh^g" "$1"

    gsed -i "s^manifold_dimension^get_manifold_dimension^g" "$1"
    gsed -i "s^space_dimension^get_space_dimension^g" "$1"
    
}
