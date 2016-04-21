import amanzi_xml.common.parameter as parameter

def valid_mesh_entity(entity):
    if entity.lower() not in ["cell", "face", "boundary_face", "node", "edge"]:
        raise RuntimeError('Unknown entity: "%s"'%entity)
    return entity.lower()
    
