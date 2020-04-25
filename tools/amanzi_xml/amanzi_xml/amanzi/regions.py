import os

from amanzi_xml import AMANZI_SRC_DIR

from amanzi_xml.utils.io import extractDoxygenXML
import amanzi_xml.utils.search as search
import amanzi_xml.common.parameter as parameter
import amanzi_xml.common.parameter_list as parameter_list
from . import mesh_entity


def createRegionLabeledSet(setname, entity, label, mesh, format="Exodus II"):
    """Create a labeled set region.

    setname     | string, name of the region
    entity      | string, entity (see mesh_entity.py)
    label       | string, label id in mesh file (note this is usually a string containing an integer)
    mesh        | string, mesh filename
    fomat       | string, format of mesh (currently only 'Exodus II' supported)

    returns the region xml
    """
    e = extractDoxygenXML(os.path.join(AMANZI_SRC_DIR, 'src', 'geometry', 'RegionLabeledSet.hh'))
    search.replace_by_name(e, "label", label)
    search.replace_by_name(e, "entity", mesh_entity.valid_mesh_entity(entity))
    search.replace_by_name(e, "format", format)
    search.replace_by_name(e, "mesh", mesh)
    pl = parameter_list.ParameterList(setname)
    pl.append(e)
    return pl

        

