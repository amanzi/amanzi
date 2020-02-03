from . import errors


#
# Functions for finding and find & replace on XML
# ---------------------------------------------------------------------------
def findall_attr(xml, attr, value):
    """Find all matches, recursively."""
    return xml.findall('.//*/[@%s="%s"]'%(attr,value))

def find_attr(xml, attr, value):
    """Find a match, recursively"""
    findall = findall_attr(xml, attr, value)
    assert len(findall) == 1
    return findall[0]

def findall_name(xml,name):
    """Find all xml objects with the given 'name'"""
    return findall_attr(xml, "name", name)
def find_name(xml,name):
    """Find an xml objects with the given 'name'"""
    return find_attr(xml, "name", name)

def findall_value(xml,value):
    """Find all xml objects with the given 'value'"""
    return findall_attr(xml, "value", value)
def find_value(xml,value):
    """Find an xml objects with the given 'value'"""
    return find_attr(xml, "value", value)

def findall_name_value(xml,name,value):
    """Find all xml objects with the given 'name' and 'value'"""
    return xml.findall('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')
def find_name_value(xml,name,value):
    """Find an xml object with the given 'name' and 'value'"""
    findall = findall_name_value(xml, name, value)
    assert len(findall) == 1
    return findall[0]

def gen_by_path(xml, names):
    assert len(names) > 0
    if len(names) == 1:
        for m in findall_name(xml,names[0]):
            yield m

    else:
        for m in findall_name(xml, names[0]):
            for n in gen_by_path(m, names[1:]): yield n

def find_by_path(xml,names):
    """Find an xml object with name path defined by list of name strings in decending hierarchical order.
    
    """
    matches = [m for m in gen_by_path(xml, names)]
    assert len(matches) == 1, "Path has %d matches"%len(matches)
    return matches[0]
  
def replace_by_value(xml,oldvalue,newvalue):
    """Replace all matches of a given 'value'='oldvalue' with 'newvalue'"""
    for r in findall_value(xml, oldvalue):
        r.setValue(newvalue)


def replace_by_name(xml,name,value):
    """Replace all matches of a given 'name' with 'newvalue'"""
    for r in findall_name(xml, name):
        r.setValue(value)

def replace_by_path(xml,names,value):
    """Replace value at end of path defined by list of name strings in decending hierarchical order."""
    find_by_path(xml,names).set('value',str(value))

def depth (xml):
    """Return depth of xml object"""
    return max([0] + [depth(child) + 1 for child in xml])
 
#
# parent map functionality
# ---------------------------------------------------------------------------
def create_parent_map(xml):
    """Creates the parent map dictionary, a map from child to parent in the XML hierarchy."""
    return {c:p for p in xml.iter() for c in p}

def replace_elem(xml, elem_sink, elem_src):
    """Replace the element 'sink' with the element 'src' in the hierarchy 'xml'"""    
    pm = create_parent_map(xml)
    p_elem_sink = pm[elem_sink]
    i = next((i for i in range(len(p_elem_sink)) if p_elem_sink[i] == elem_sink))
    p_elem_sink[i] = elem_src

def get_parent(xml,elem,level=1):
    """Parses up the parent map by 'level'"""
    pm = create_parent_map(xml)
    parent = elem
    for i in range(level):
        parent = pm[parent]
    return parent

def get_path(xml, elem, level=None):
    """Parses up the parent map through the entire hierarchy.

    Returns a list, [xml, ..., elem_parent, elem]
    """
    pm = create_parent_map(xml)
    parent = elem
    path = []
    path.append(parent)
    if level is None:
        while parent != xml:
            parent = pm[parent]
            path.append(parent)
    else:
        counter = 0
        while parent != xml and counter < level:
            parent = pm[parent]
            path.append(parent)
            counter += 1
    
    path.reverse()
    return path

def get_path_namelist(xml, elem, level=None):
    """Parses up the parent map through the entire hierarchy and returns list of path names.

    Returns a list, [xml, ..., elem_parent, elem]
    """
    return [e.attrib['name'] for e in get_path(xml,elem,level)]

def print_path(xml,elem,level=None):
    """Parses up the parent map through the entire hierarchy and prints to terminal.
    """
    elems = get_path(xml,elem,level)
    ind = ''
    s = ''
    for e in elems:
        s += "\n"+ind+e.attrib['name']
        ind += '    '
    print(s)

def get_value(xml, name):
    """Return value associated with name

    Returns a single value if there is one occurrence of name,
    a list if there is more than one.
    """
    out = []
    for r in findall_name(xml,name):
        out.append(r.attrib['value'])
    if len(out) > 1: return out
    else: return out[0]

def global_remove(xml, elem):
    """Removes the xml 'elem' wherever it is locaed in 'xml'"""
    pm = create_parent_map(xml)
    pm[elem].remove(elem)




#
# FIX AND REMOVE ME!
# --------------------
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

