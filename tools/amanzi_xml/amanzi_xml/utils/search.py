import io

# Searches based upon the "name" attribute as the unique identifyer.
def _getChildByName(elem, name):
    """Assumed that names are unique!"""
    for subel in elem:
        if (name == subel.get("name")):
            return subel

# Searches based upon tags -- no uniqueness!
def _tagGenerator(elem, tag):
    for subel in elem:
        if tag == subel.tag:
            yield subel

def getElementByName(elem, path):
    """Assumed that names are unique!"""
    enames = path.strip("/").split("/")
    assert(elem.get("name") == enames[0])
    enames.pop(0)

    if len(enames) == 0:
        return elem
    else:
        return getElementByName(_getChildByName(elem, enames[0]), "/".join(enames))

def getElementByPath(elem, path):
    if elem.tag == "ParameterList":
        return getElementByName(elem, path)

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
        for child in _tagGenerator(elem, childtagname[0]):
            if len(childtagname) == 1:
                return getElementByPath(child, "/".join(etagnames))
            else:
                if child.get("name") == childtagname[1]:
                    return getElementByPath(child, "/".join(etagnames))


# def searchAndRemoveByName(pl, abspath):
#     """Search for an absolute path and remove the parameter."""
#     subpath = abspath.split("/")
#     containing_path = "/".join(subpath[:-1])
#     container = getElementByName(pl, containing_path)
#     container.remove(container.getElement(subpath[-1]))

# def searchAndReplaceByName(pl, changeset):
#     """Search for a path and replace the value.

#     Changeset is expected of the form:
#       path/to/my/parameter=newvalue

#     or
#       /abs/path/to/my/parameter=newvalue
#     """
#     split = changeset.split("=")
#     if len(split) != 2:
#         raise RuntimeError("Invalid changeset %s not of form path=val"%changeset)

#     path,val = tuple(split)
#     if path.startswith("/"):
#         elem = getElementByName(pl, path)
#         elem.setValue(val)
#     else:
#         for elem in findElementNameGenerator(pl,path):
#             elem.setValue(val)

