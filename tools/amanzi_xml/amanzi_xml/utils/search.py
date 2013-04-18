import io

def _findSingleElementGenerator(elem, name):
    if elem.get("name") == name:
        yield elem
    else:
        for subel in elem:
            for entry in _findSingleElementGenerator(subel, name):
                yield entry

def findElementGenerator(elem, path):
    """Generator that parses through all Elements with given relative path/name.

    The path argument may be either a name or a path, which may be
    either absolute or relative, i.e.:
      findElementGenerator(elem, "Mesh")
      findElementGenerator(elem, "Mesh/Expert/Verify Mesh")
      findElementGenerator(elem, "/Mesh/Expert/Verify Mesh")

    """
    enames = path.split("/")
    if len(enames) == 1:
        for entry in _findSingleElementGenerator(elem, enames[0]):
            yield entry
    else:
        for entry in _findSingleElementGenerator(elem, enames[0]):
            for subentry in findElementGenerator(entry, "/".join(enames[1:])):
                yield subentry

def getElementByPath(elem, path):
    """Generator that parses through all Elements from an absolute path"""
    if not path.startswith("/"):
        raise ValueError("getElementByPath must be provided an absolute path")

    enames = path.strip("/").split("/")
    cur_elem = elem
    first = enames.pop(0)
    assert cur_elem.get("name") == first

    while len(enames) > 0:
        nextname = enames.pop(0)
        cur_elem = cur_elem.getElement(nextname)

    return cur_elem


def searchAndReplace(pl, changeset):
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
        elem = getElementByPath(pl, path)
        elem.setValue(val)
    else:
        for elem in findElementGenerator(pl,path):
            elem.setValue(val)

def searchAndRemove(pl, abspath):
    """Search for an absolute path and remove the parameter."""
    subpath = abspath.split("/")
    containing_path = "/".join(subpath[:-1])
    container = getElementByPath(pl, containing_path)
    container.remove(container.getElement(subpath[-1]))
