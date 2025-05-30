import xml.etree.ElementTree as etree
from . import errors as aerrors
from . import io as aio

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
            print((elem, cname))
            print((list(objects.keys())))
            raise aerrors.NotNativeSpecError()
        else:
            return obj.from_Element(elem)

    cname = elem.get('name')
    if cname in objects:
        # a typed list
        return objects[cname].from_Element(elem)
    else:
        # a genericl parameter list
        return objects["ParameterList"].from_Element(elem)


def fromElement(elem):
    """Creates a amanzi-xml hierarchy from standard Element"""
    if elem.tag == 'Parameter':
        return _parameterFromElement(elem)
    elif elem.tag == 'ParameterList':
        return _listObjectFromElement(elem)
    elif 'function Comment' in str(elem.tag):
        return objects['function Comment'].from_Element(elem)
    raise RuntimeError('Invalid element with tag %s' % elem.tag)


