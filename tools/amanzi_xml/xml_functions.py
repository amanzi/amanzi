"""Collection of useful utility functions for working with XML objects"""

try:
    import xml.etree.cElementTree as ET
except:
    import xml.etree.ElementTree as ET

#
# Functions for adding sub-entries to the X
# ---------------------------------------------------------------------------
def add_Parameter(xml, name, value, vtype):
    """Add a Parameter to the input xml.

    Little error checking here.
    """
    if type(value) is list:
        value = "{"+",".join(str(v) for v in value)+"}"
    ET.SubElement(xml,'Parameter',{'name':name,'type':vtype,'value':str(value)})

def add_ParameterList(xml, name):
    """Add a ParameterList to the input xml.

    Little error checking here.
    """
    ET.SubElement(xml,'ParameterList',{'name':name})
    

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

def replace_by_value(xml,oldvalue,newvalue):
    """Replace all matches of a given 'value'='oldvalue' with 'newvalue'"""
    for r in findall_value(xml, oldvalue):
        r.set('value', str(newvalue))

def replace_by_name(xml,name,value):
    """Replace all matches of a given 'name' with 'newvalue'"""
    for r in findall_name(xml, name):
        r.set('value', str(value))

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
    i = (i for i in range(len(p_elem_sink)) if p_elem_sink[i] == elem_sink).next()
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

def remove(xml, elem):
    """Removes the xml 'elem' wherever it is locaed in 'xml'"""
    pm = create_parent_map(xml)
    pm[elem].remove(elem)
