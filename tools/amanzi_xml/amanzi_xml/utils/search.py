import io
import errors

# Searches based upon the "name" attribute as the unique identifier.
def generateChildByName(elem, name):
    """Generator for children of a given name"""
    #print "generateChildByName(%s, %s)"%(elem.get("name"), name)
    for subel in elem:
        if (name == subel.get("name")):
            yield subel

def generateChildByNamePath(elem, path):
    """Generator for (grand)children of a given path"""
    #print "generateChildByNamePath(%s, %s)"%(elem.get("name"), path)
    enames = path.strip("/").split("/")
    assert len(enames) > 0
    if len(enames) == 1:
        for match in generateChildByName(elem, enames[0]):
            yield match
    else:
        for match in generateChildByName(elem, enames[0]):
            for submatch in generateChildByNamePath(match, "/".join(enames[1:])):
                yield submatch


def childByName(elem, name):
    """Gets child of a given name"""
    children = [el for el in generateChildByName(elem, name)]
    if len(children) == 0:
        raise errors.MissingXMLError()
    elif len(children) > 1:
        raise errors.NonUniqueXMLError()
    else:
        return children[0]

def childByNamePath(elem, path):
    """Gets a (grand)child of a given name"""
    children = [el for el in generateChildByNamePath(elem, path)]
    if len(children) == 0:
        raise errors.MissingXMLError()
    elif len(children) > 1:
        raise errors.NonUniqueXMLError()
    else:
        return children[0]

def getElementByNamePath(elem, path):
    """Wrapper for childByNamePath which checks the current elem"""
    enames = path.strip("/").split("/")
    assert len(enames) > 0
    if elem.get("name") == enames[0]:
        enames.pop(0)
    if len(enames) == 0:
        return elem
    else:
        return childByNamePath(elem, "/".join(enames))

def generateElementByNamePath(elem, path):
    """Searches for all (grand)children that match the relative namepath"""
    #print "generateElementByNamePath(%s, %s)"%(elem.get("name"), path)
    enames = path.strip("/").split("/")
    assert len(enames) > 0

    def generateSingleElementByName(elem, name):
        #print "generateSingleElementByName(%s, %s)"%(elem.get("name"), name)
        if elem.get("name") == name:
            yield elem
        else:
            for el in elem:
                for match in generateSingleElementByName(el, name):
                    yield match

    if len(enames) == 1:
        for match in generateSingleElementByName(elem, enames[0]):
            yield match
    else:
        for match in generateSingleElementByName(elem, enames[0]):
            for submatch in generateChildByNamePath(match, "/".join(enames[1:])):
                yield submatch


# Searches based upon tags -- no uniqueness!
def generateChildByTag(elem, tag):
    for subel in elem:
        if tag == subel.tag:
            yield subel


def getElementByTagPath(elem, path):
    """Get element by tag, which is the amanzi_xml-aware ParameterList used for creating specs"""
    if elem.tag == "ParameterList":
        return getElementByNamePath(elem, path)

    etagnames = path.strip("/").split("/")
    mytagname = etagnames[0].strip().strip("}").strip("{").split(",")

    if len(mytagname) > 0:
        assert elem.tag == mytagname[0]
    if len(mytagname) > 1:
        assert elem.get("name") == mytagname[1]
    etagnames.pop(0)

    if len(etagnames) == 0:
        return elem
    else:
        childtagname = etagnames[0].strip().strip("}").strip("{").split(",")
        for child in generateChildByTag(elem, childtagname[0]):
            if len(childtagname) == 1:
                return getElementByNamePath(child, "/".join(etagnames))
            else:
                if child.get("name") == childtagname[1]:
                    return getElementByNamePath(child, "/".join(etagnames))


def getElementByTags(elem, path):
    """Get element by tags only, which is for v2 XML specs"""
    etagnames = path.strip("/").split("/")
    mytagname = etagnames[0].strip().strip("}").strip("{").split(",")

    if len(mytagname) > 0:
        assert elem.tag == mytagname[0]
    etagnames.pop(0)

    if len(etagnames) == 0:
        return elem
    else:
        childtagname = etagnames[0].strip().strip("}").strip("{").split(",")
        for child in generateChildByTag(elem, childtagname[0]):
            if len(childtagname) == 1:
                return getElementByTags(child, "/".join(etagnames))
            else:
                if child.get("name") == childtagname[1]:
                    return getElementByTags(child, "/".join(etagnames))


# def searchAndRemoveByName(pl, abspath):
#     """Search for an absolute path and remove the parameter."""
#     subpath = abspath.split("/")
#     containing_path = "/".join(subpath[:-1])
#     container = getElementByName(pl, containing_path)
#     container.remove(container.getElement(subpath[-1]))

def searchAndReplaceByName(pl, changeset):
    """Search for a path and replace the value.

    Changeset is expected of the form:
      path/to/my/parameter=newvalue

    or
      /abs/path/to/my/parameter=newvalue
    """
    split = changeset.split("=")
    if len(split) != 2:
        raise RuntimeError("Invalid changeset %s not of form path=val"%changeset)

    path,val = tuple(split)
    if path.startswith("/"):
        elem = getElementByNamePath(pl, path)
        elem.setValue(val)
    else:
        for elem in generateElementByNamePath(pl,path):
            elem.setValue(val)

