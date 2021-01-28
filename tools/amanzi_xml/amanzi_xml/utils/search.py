""" Includes functions for searching xml objects by name and path."""

from . import errors
import warnings

#
# Functions for finding elements based on their own attributes
# ---------------------------------------------------------------------------
def children_by_tag(xml, tag):
    """Find all matches, at depth 1.

    Note that tag may be a tag, or it may be "tag,name".  This is ugly
    but useful.
    """
    name = None
    if ',' in tag:
        name = ','.join(tag.split(',')[1:])
        tag = tag.split(',')[0]

    matches = xml.findall(tag)
    if name is not None:
        matches = [m for m in matches if m.get("name") == name]
    return matches
        
def child_by_tag(xml, tag):
    """Find unique child by tag.  Note tags may not be unique."""
    matches = children_by_tag(xml, tag)
    if len(matches) == 0:
        raise errors.MissingXMLError(f"Object does not have child with tag '{tag}'")
    elif len(matches) > 1:
        raise errors.NonUniqueXMLError(f"Object has more than one child with tag '{tag}'")
    return matches[0]

def children_by_attr(xml, attr, value):
    """Find all matches, at depth 1."""
    return xml.findall(f'./*[@{attr}="{value}"]')

def children_by_name(xml, name):
    """Find children by name."""
    return children_by_attr(xml, 'name', name)

def child_by_name(xml, name):
    """Find unique child by name.  Note names should always be unique."""
    matches = children_by_name(xml, name)
    if len(matches) == 0:
        raise errors.MissingXMLError(f"Object does not have child '{name}'")
    elif len(matches) > 1:
        raise errors.NonUniqueXMLError(f"Object has more than one child '{name}'")
    return matches[0]

# ---------------------------------------------------------------------------


def findall_attr(xml, attr, value):
    """Find all matches, at any depth, in the xml tree."""
    return xml.findall('.//*/[@{}="{}"]'.format(attr,value))

def find_attr(xml, attr, value):
    """Find a unique xml object, at any depth, in the xml tree."""
    matches = findall_attr(xml, attr, value)
    if len(matches) == 0:
        raise errors.MissingXMLError(f"Object does not have grandchild '{attr}={value}'")
    elif len(matches) > 1:
        raise errors.NonUniqueXMLError(f"Object has more than one grandchild '{attr}={value}'")
    return matches[0]

def findall_name(xml, name):
    """Find all xml objects with the given 'name' at any depth."""
    return findall_attr(xml, "name", name)

def find_name(xml, name):
    """Find a unique xml object with the given 'name' at any depth."""
    return find_attr(xml, "name", name)

def findall_value(xml, value):
    """Find all xml objects with the given 'value'"""
    return findall_attr(xml, "value", value)

def find_value(xml,value):
    """Find an xml objects with the given 'value'"""
    return find_attr(xml, "value", value)

def findall_name_value(xml, name, value):
    """Find all xml objects at any depth with the given 'name' and 'value'"""
    return xml.findall('.//*/[@name="'+str(name)+'"]/[@value="'+str(value)+'"]')

def find_name_value(xml, name, value):
    """Find a unique xml object at any depth with the given 'name' and 'value'"""
    matches = findall_name_value(xml, name, value)
    if len(matches) == 0:
        raise errors.MissingXMLError(f"Object does not have grandchild '{name},{value}'")
    elif len(matches) > 1:
        raise errors.NonUniqueXMLError(f"Object has more than one grandchild '{name},{value}'")
    return matches[0]

#
# Finds by path, a list of tags or names for nesting
# ---------------------------------------------------------------------------
def gen_by_tag_path(xml, tags):
    """Generator that takes a list of tags and returns matching elements."""
    assert(type(tags) is not str)
    assert(len(tags) > 0)
    if len(tags) == 1:
        for m in children_by_tag(xml,tags[0]):
            yield m
    else:
        for m in children_by_tag(xml,tags[0]):
            for n in gen_by_tag_path(m, tags[1:]):
                yield n


def gen_by_path(xml, names, no_skip=False):
    """Generator that takes a list of names and returns matching elements.

    If no_skip is True, successive names must be direct children.
    """
    if no_skip:
        findall = children_by_name
    else:
        findall = findall_name

    assert(type(names) is not str)
    assert(len(names) > 0)
    if len(names) == 1:
        for m in findall(xml,names[0]):
            yield m
    else:
        for m in findall(xml, names[0]):
            for n in gen_by_path(m, names[1:], no_skip):
                yield n

def find_tag_path(xml, tags):
    """Find a unique xml object with a list of tags."""
    assert(type(tags) is not str)
    assert(len(tags) > 0)
    if (xml.tag == tags[0]):
        tags.pop(0)

    matches = list(gen_by_tag_path(xml, tags))
    if len(matches) == 0:
        raise errors.MissingXMLError(f"Object does not have tag path '{tags}'")
    elif len(matches) > 1:
        raise errors.NonUniqueXMLError(f"Object has more than one tag path '{tags}'")
    return matches[0]

def findall_tag_path(xml, tags):
    """Find a list of matching elements by tag."""
    assert(type(tags) is not str)
    assert(len(tags) > 0)
    if (xml.tag == tags[0]):
        tags.pop(0)
    return list(gen_by_tag_path(xml, tags))

def find_path(xml, names, no_skip=False):
    """Find a unique xml object with a list of names."""
    matches = list(gen_by_path(xml, names, no_skip))
    if len(matches) == 0:
        raise errors.MissingXMLError("Object does not have path '{}'".format(names))
    elif len(matches) > 1:
        raise errors.NonUniqueXMLError("Object has more than one path '{}'".format(names))
    return matches[0]

def findall_path(xml, names, no_skip=False):
    """Find a list of matching elements."""
    return list(gen_by_path(xml, names, no_skip))

    
#
# parent map functionality
# ---------------------------------------------------------------------------
def global_depth(xml):
    """Return maximal depth of an xml tree."""
    return max([0] + [global_depth(child) + 1 for child in xml])

def depth(xml1, xml2):
    """Assuming xml2 is a child of xml1, calculates how deep."""
    pm = parent_map(xml1)
    d = 0
    while xml2 is not xml1:
        d += 1
        xml2 = pm[xml2]
    return d
 
def parent_map(xml):
    """Creates the parent map dictionary, a map from child to parent in the XML hierarchy."""
    try:
        return xml.pm
    except AttributeError:
        xml.pm = {c:p for p in xml.iter() for c in p}
    return xml.pm

def get_parent(xml,elem,level=1):
    """Parses up the parent map by 'level'"""
    pm = parent_map(xml)
    parent = elem
    for i in range(level):
        parent = pm[parent]
    return parent

def get_path(xml, elem, level=None):
    """Parses up the parent map through the entire hierarchy.

    Returns a list, [xml, ..., elem_parent, elem]
    """
    pm = parent_map(xml)
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
    
    return list(reversed(path))

def get_path_namelist(xml, elem, level=None):
    """Parses up the parent map through the entire hierarchy and returns list of path names.

    Returns a list, [xml, ..., elem_parent, elem]
    """
    return [e.get('name') for e in get_path(xml, elem, level)]

    
#
# high level search and replace
#
def change_value(xml, name_or_names, value, allow_multiple=False, no_skip=False):
    """Changes the value of an element in xml.

    Parameters
    ----------
    xml : etree.Element
      The parent list to search.
    name_or_names : str or list(str)
      Either a single name or list of names defining a path to the object to
      change.
    value : str,int,float,list,etc
      The new value.
    allow_multiple : bool, optional=False
      Allow and change multiple matches.
    no_skip : bool, optional=False
      If True, successive names must be direct children.
    """
    if type(name_or_names) is str:
        name_or_names = [name_or_names,]

    matches = findall_path(xml, name_or_names, no_skip)
    if len(matches) == 0:
        raise errors.MissingXMLError("Object does not have path '{}'".format(name_or_names))

    if not allow_multiple and len(matches) > 1:
        raise errors.NonUniqueXMLError("Object has more than one path '{}'".format(name_or_names))

    for match in matches:
        match.setValue(value)


def change_name(xml, old_name_or_names, new_name, allow_multiple=False, no_skip=False):
    """Changes the name of an element in xml.

    Parameters
    ----------
    xml : etree.Element
      The parent list to search.
    old_name_or_names : str or list(str)
      Either a single name or list of names defining a path to the elemnt(s) to
      change.
    new_name : str
      The new name of the found element(s).
    allow_multiple : bool, optional=False
      Allow and change multiple matches.
    no_skip : bool, optional=False
      If True, successive names must be direct children.
    """
    if type(old_name_or_names) is str:
        old_name_or_names = [old_name_or_names,]

    matches = findall_path(xml, old_name_or_names, no_skip)
    if len(matches) == 0:
        raise errors.MissingXMLError("Object does not have path '{}'".format(old_name_or_names))

    if not allow_multiple and len(matches) > 1:
        raise errors.NonUniqueXMLError("Object has more than one path '{}'".format(old_name_or_names))

    for match in matches:
        match.setName(new_name)
        

def remove_element(xml, name_or_names, allow_multiple=False, no_skip=False):
    """Removes an element by name or path.

    Parameters
    ----------
    xml : etree.Element
      The parent list to search.
    name_or_names : str or list(str)
      Either a single name or list of names defining a path to the object to
      change.
    allow_multiple : bool, optional=False
      Allow and change multiple matches.
    no_skip : bool, optional=False
      If True, successive names must be direct children.

    Returns
    -------
    old_el : the removed element(s)
    """
    if type(name_or_names) is str:
        name_or_names = [name_or_names,]

    matches = findall_path(xml, name_or_names, no_skip)
    if len(matches) == 0:
        raise errors.MissingXMLError("Object does not have path '{}'".format(name_or_names))

    if not allow_multiple and len(matches) > 1:
        raise errors.NonUniqueXMLError("Object has more than one path '{}'".format(names))

    pm = parent_map(xml)
    removed = []
    for match in matches:
        parent = pm[match]
        removed.append(match)
        parent.remove(match)
    if allow_multiple:
        return removed
    else:
        assert(len(removed) == 1)
        return removed[0]

def replace_element(xml, name_or_names, new_el_or_els, allow_multiple=False, no_skip=False):
    """Replaces an element with another 

    Parameters
    ----------
    xml : etree.Element
      The parent list to search.
    name_or_names : str or list(str)
      Either a single name or list of names defining a path to the object to
      change.
    new_el_or_els : etree.Element or list(etree.Element)
      The new element(s) to add.
    allow_multiple : bool, optional=False
      Allow and change multiple matches.
    no_skip : bool, optional=False
      If True, successive names must be direct children.

    Returns
    -------
    old_el : the removed element(s)
    """
    if type(name_or_names) is str:
        name_or_names = [name_or_names,]
    if type(new_el_or_els) is not list:
        new_el_or_els = [new_el_or_els,]

    matches = findall_path(xml, name_or_names, no_skip)
    if len(matches) == 0:
        raise errors.MissingXMLError("Object does not have path '{}'".format(names))

    if not allow_multiple and len(matches) > 1:
        raise errors.NonUniqueXMLError("Object has more than one path '{}'".format(names))

    pm = parent_map(xml)
    removed = []
    for match in matches:
        parent = pm[match]
        removed.append(parent.pop(match))
        parent.extend(new_el_or_els)
    if allow_multiple:
        return removed
    else:
        assert(len(removed) == 1)
        return removed[0]
        
    
    
