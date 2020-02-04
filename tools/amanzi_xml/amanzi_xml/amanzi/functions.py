import os

from amanzi_xml import AMANZI_SRC_DIR

from amanzi_xml.utils.io import extractDoxygenXML
import amanzi_xml.utils.search as search
import amanzi_xml.common.parameter as parameter
import amanzi_xml.common.parameter_list as parameter_list
from . import mesh_entity


def createConstantFunction(value):
    """Create a constant function with with a given value, region, and entity.

    value       | double, value of the function
    """
    e = extractDoxygenXML(os.path.join(AMANZI_SRC_DIR, 'src', 'functions', 'ConstantFunction.hh'))
    search.replace_by_name(e, "value", value)
    return e

def createFunctionSpec(spec_name, region, component, func):
    """Create a Function spec on a 'region' and 'component'"""
    if type(region) is str:
        region_p = parameter.StringParameter("region", region)
    else:
        region_p = parameter.ArrayStringParameter("regions", region)

    if type(component) is str:
        component_p = parameter.StringParameter("component", component)
    else:
        component_p = parameter.ArrayStringParameter("components", component)

    pl = parameter_list.ParameterList(spec_name)
    pl.append(region_p)
    pl.append(entity_p) # TODO: entity_p undefined: did you mean component_p?
    pl.sublist("function").append(func)
    return pl
        

