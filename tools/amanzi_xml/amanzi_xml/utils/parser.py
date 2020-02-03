import xml.etree.ElementTree as etree
from . import errors as aerrors
import io as aio

objects = dict()

class ValidationLevel(dict):
    def __init__(self, missing_required=True,
                       type_required=True,
                       optional_type=True,
                       existing_not_optional=False):
        self['missing_required'] = missing_required
        self['type_required'] = type_required
        self['optional_type'] = optional_type
        self['existing_not_optional'] = existing_not_optional
        if missing_required is None and type_required is None and optional_type is None and existing_not_optional is None:
            self.none = True
        else:
            self.none = False



def _parameterFromElement(elem):
    return objects[elem.get("type")].from_Element(elem)

def _listObjectFromElement(elem):
    # see if this is a new-style elem, whose Lists have Types
    cname = elem.get("type")
    if cname is not None:
        try:
            obj = objects[cname]
        except KeyError:
            print("Error interpreting list:")
            print(elem, cname)
            print(objects.keys())
            raise aerrors.NotNativeSpecError()
        else:
            return obj.from_Element(elem)

    # else try to guess the type from the List name
    cname = elem.get('name')
    assert cname is not None
    try:
        return objects[cname].from_Element(elem)

    # else just create a generic ParameterList
    except KeyError:
        return objects["ParameterList"].from_Element(elem)


def fromElement(elem):
    """Creates a amanzi-xml hierarchy from standard Element"""
    if elem.tag == 'Parameter':
        return _parameterFromElement(elem)
    elif elem.tag == 'ParameterList':
        return _listObjectFromElement(elem)
    else:
        raise RuntimeError('Invalid element with tag %s'%elem.tag)


