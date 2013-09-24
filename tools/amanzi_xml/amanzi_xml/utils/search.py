import io

# Searches based upon the "name" attribute as the unique identifyer.
def _findSingleElementNameGenerator(elem, name):
    if elem.get("name") == name:
        yield elem
    else:
        for subel in elem:
            for entry in _findSingleElementNameGenerator(subel, name):
                yield entry

def findElementNameGenerator(elem, path):
    """Generator that parses through all Elements with given relative path/name.

    The path argument may be either a name or a path, which may be
    either absolute or relative, i.e.:
      findElementNameGenerator(elem, "Mesh")
      findElementNameGenerator(elem, "Mesh/Expert/Verify Mesh")
      findElementNameGenerator(elem, "/Mesh/Expert/Verify Mesh")

    """
    enames = path.split("/")
    if len(enames) == 1:
        for entry in _findSingleElementNameGenerator(elem, enames[0]):
            yield entry
    else:
        for entry in _findSingleElementNameGenerator(elem, enames[0]):
            for subentry in findElementNameGenerator(entry, "/".join(enames[1:])):
                yield subentry

def getElementByName(elem, path):
    """Generator that parses through all Elements from an absolute path"""
    if not path.startswith("/"):
        raise ValueError("getElementByName must be provided an absolute path")

    enames = path.strip("/").split("/")
    cur_elem = elem
    first = enames.pop(0)
    assert cur_elem.get("name") == first

    while len(enames) > 0:
        nextname = enames.pop(0)
        cur_elem = cur_elem.getElement(nextname)

    return cur_elem


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
        elem = getElementByName(pl, path)
        elem.setValue(val)
    else:
        for elem in findElementNameGenerator(pl,path):
            elem.setValue(val)

def searchAndRemoveByName(pl, abspath):
    """Search for an absolute path and remove the parameter."""
    subpath = abspath.split("/")
    containing_path = "/".join(subpath[:-1])
    container = getElementByName(pl, containing_path)
    container.remove(container.getElement(subpath[-1]))



# Searches based upon tags -- no uniqueness!
def _findSingleElementTagGenerator(elem, tag):
    if elem.tag == tag:
        yield elem
    else:
        for subel in elem:
            for entry in _findSingleElementTagGenerator(subel, tag):
                yield entry

def findElementTagGenerator(elem, path):
    """Generator that parses through all Elements with given relative path/tag.

    The path argument may be either a tag or a path, which may be
    either absolute or relative, i.e.:
      findElementTagGenerator(elem, "Mesh")
      findElementTagGenerator(elem, "Mesh/Expert/Verify Mesh")
      findElementTagGenerator(elem, "/Mesh/Expert/Verify Mesh")

    """
    etags = path.split("/")
    if len(etags) == 1:
        for entry in _findSingleElementTagGenerator(elem, etags[0]):
            yield entry
    else:
        for entry in _findSingleElementTagGenerator(elem, etags[0]):
            for subentry in findElementTagGenerator(entry, "/".join(etags[1:])):
                yield subentry


# Common interface
def getElementByPath(xml, path):
    if xml.tag == "ParameterList" or xml.tag == "Parameter":
        return findElementNameGenerator(xml, path).next()
    else:
        return findElementTagGenerator(xml, path).next()
