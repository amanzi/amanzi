################################################################################

import sys, types, os

from xml.etree.ElementTree import Comment, ProcessingInstruction, QName
from xml.etree.ElementTree import _encode, _escape_cdata, _escape_attrib
from xml.etree.ElementTree import ElementTree, Element, _ElementInterface

################################################################################
'''
Global list arrays for the different data types
'''
_bool_type_strings = ['bool', 'boolean']
_int_type_strings  = ['int', 'integer']
_float_type_strings = ['real', 'double', 'float']
_string_type_strings = ['str', 'string', 'char', 'character']

################################################################################

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

################################################################################

def prettyprint(tree,file=sys.stdout,encoding='utf-8',sortflag='default',sortcmp=None):
    root = tree.getroot()
    indent(root)
    tree.write(file=file,encoding=encoding,sortflag=sortflag,sortcmp=sortcmp)

################################################################################

def getTypeFromString(str_type):
    test_str = str_type.lower
    if test_str in _bool_type_strings:
        return types.BooleanType
    elif test_str in _int_type_strings:
        return types.IntType
    elif test_str in _float_type_strings:
        return types.FloatType
    elif test_str in _string_type_strings:
        return types.StringType

################################################################################

def getStringFromArray(a,delim=','):
    return delim.join([str(s) for s in a])

################################################################################

def getArrayFromString(buf,type='string',delim=','):
    test = type.lower()
    if test in _int_type_strings :
        return [int(s) for s in buf.split(delim)]
    elif test in _float_type_strings :
        return [float(s) for s in buf.split(delim)]
    else:
        return [str(s) for s in buf.split(delim)]
    
################################################################################

class InputTree(ElementTree):

    # Code Name
    code_name = None

    # Root Element
    root = None

    def __init__(self,root_tag='AmanziInput',code_name=None):
        ElementTree.__init__(self)

        # Define the root element
        root = Element(root_tag)
        if code_name is not None:
            root.set('code_name',code_name)

        self._setroot(root)
        self.root = root


    def attach(self,node):
        self.root.append(node)

    def dumpXML(self,file=sys.stdout,encoding='utf-8',xml_translate=None,*args):

        sortflag='name,label,type'
        if xml_translate is not None :
            print_tree = xml_translate(self,*args)
            prettyprint(print_tree,file,encoding,sortflag=sortflag)
        else:
            sortflag = sortflag + ',vector,length,delim'
            prettyprint(self,file,encoding,sortflag=sortflag)

            ##
    # Writes the element tree to a file, as XML.
    #
    # @param file A file name, or a file object opened for writing.
    # @param encoding Optional output encoding (default is US-ASCII).

    def write(self, file, encoding="us-ascii", sortflag="default", sortcmp=None):
        assert self._root is not None
        if not hasattr(file, "write"):
            file = open(file, "wb")
        if not encoding:
            encoding = "us-ascii"
        elif encoding != "utf-8" and encoding != "us-ascii":
            file.write("<?xml version='1.0' encoding='%s'?>\n" % encoding)
        self._write(file, self._root, encoding, {}, sortflag, sortcmp)

    def _write(self, file, node, encoding, namespaces, sortflag="default", sortcmp=None): # don't break existing code that relies on _write()s parameters, if any
        # write XML to file
        tag = node.tag
        if tag is Comment:
            file.write("<!-- %s -->" % _escape_cdata(node.text, encoding))
        elif tag is ProcessingInstruction:
            file.write("<?%s?>" % _escape_cdata(node.text, encoding))
        else:
            items = node.items()
            xmlns_items = [] # new namespaces in this scope
            try:
                if isinstance(tag, QName) or tag[:1] == "{":
                    tag, xmlns = fixtag(tag, namespaces)
                    if xmlns: xmlns_items.append(xmlns)
            except TypeError:
                _raise_serialization_error(tag)
            file.write("<" + _encode(tag, encoding))
            if items or xmlns_items:
                
                ##NEW

                if sortflag!="default":
                    if ":" not in sortflag:
                        sortflag = ":"+sortflag
                    sortitems = sortflag.split(";")
                    try:
                        sortitems = [[tagdef.split(":")[0], tagdef.split(":")[1].split(",")] for tagdef in sortitems]
                        temp = []
                        for tagdef in sortitems:
                            if tagdef[0] in node.tag:
                                for sortitem in tagdef[1]:
                                    temp.extend([i for i in items if sortitem in i]) # add what matches pattern
                                break
                        # then sort and add what's left
                        items.sort(cmp=sortcmp)
                        temp.extend([i for i in items if i not in temp])
                        items = temp
                    except IndexError:
                        sys.stderr.write("sortflag not formatted correctly, order won't be applied")
                        try:
                            items.sort(cmp=sortcmp)
                        except Exception:
                            sys.stderr.write("sortcmp not a valid comparator, sorting alphabetically instead")
                            items.sort()
                else:
                    try:
                        items.sort(cmp=sortcmp)
                    except Exception:
                        sys.stderr.write("sortcmp not a valid comparator, sorting alphabetically instead")
                        items.sort()
                ###

                        
                for k, v in items:
                    try:
                        if isinstance(k, QName) or k[:1] == "{":
                            k, xmlns = fixtag(k, namespaces)
                            if xmlns: xmlns_items.append(xmlns)
                    except TypeError:
                        _raise_serialization_error(k) # TODO: undefined
                    try:
                        if isinstance(v, QName):
                            v, xmlns = fixtag(v, namespaces) # TODO: undefined
                            if xmlns: xmlns_items.append(xmlns)
                    except TypeError:
                        _raise_serialization_error(v) # TODO: undefined
                    file.write(" %s=\"%s\"" % (_encode(k, encoding),
                                               _escape_attrib(v, encoding)))
                for k, v in xmlns_items:
                    file.write(" %s=\"%s\"" % (_encode(k, encoding),
                                               _escape_attrib(v, encoding)))
            if node.text or len(node):
                file.write(">")
                if node.text:
                    file.write(_escape_cdata(node.text, encoding))
                for n in node:
                    self._write(file, n, encoding, namespaces, sortflag, sortcmp)
                file.write("</" + _encode(tag, encoding) + ">")
            else:
                file.write(" />")
            for k, v in xmlns_items:
                del namespaces[v]
        if node.tail:
            file.write(_escape_cdata(node.tail, encoding))
        
################################################################################

class _InputElementInterface(_ElementInterface):

    # Name of the Input Element node
    name = None

    # Type of element either Parameter or a Branch
    elemtype = None

    # 
    def __init__(self,name,elemtype):
        attr = {}
        attr['name'] = name
        _ElementInterface.__init__(self,elemtype,attr)
        self.name = name
        self.elemtype = elemtype


    def getName(self):
        return self.name

    def elemType(self):
        return self.elemtype

    def isaParameterElement(self):
        if self.elemType == 'Parameter':
            return True
        else:
            return False

    def isBranchElement(self):
        if self.elemType == 'Branch':
            return True
        else:
            return False

################################################################################

class _ParameterElementInterface(_InputElementInterface):

    # Python Data Type 
    py_type = types.NoneType

    # Dictionary that maps string type value to a Python data type
    _type_map = {}

    def __init__(self,name,type,data,delim=','):
        _InputElementInterface.__init__(self,name,'Parameter')
       
        '''
        Set the type, two kinds, one the user defined type and
        the internal Python built-in type
        '''
        for s in _bool_type_strings:
            self._type_map[s] = types.BooleanType
        for s in _int_type_strings:
            self._type_map[s] = types.IntType
        for s in _float_type_strings:
            self._type_map[s] = types.FloatType
        for s in _string_type_strings:
            self._type_map[s] = types.StringType

        if type in self._type_map:
            self.py_type = self._type_map[type]
        else:
            self.py_type = types.NoneType

        self.set('type',type)

        if isinstance(data,list):
            list_len = str(len(data))
            self.set('vector','True')
            self.set('length',list_len)
            self.set('delim',delim)
            self.text = getStringFromArray(data,delim)
        else:
            self.text = str(data)

    def getPythonType(self):
        return self.py_type

    def getStringType(self):
        return self.type

################################################################################

class _BranchElementInterface(_InputElementInterface):

    def __init__(self,name):
        _InputElementInterface.__init__(self,name,'Branch')

################################################################################

def InputElement(name,elemtype):
    return _InputElementInterface(name,elemtype)

################################################################################

def ParameterElement(name,type,data,delim=','): 
    return _ParameterElementInterface(name,type,data,delim)

################################################################################

def BranchElement(name): 
    return _BranchElementInterface(name)
    


# If run as a script, do some testing
if __name__ == '__main__':

    # Set up some input parameters and dump to STDOUT
    input = InputTree()

    # Scalar Parameters
    node = ParameterElement('Title','string','<< The Title >>')
    input.attach(node)

    node = ParameterElement('Permability','double', 1.234e-5)
    input.attach(node)

    vec = [1.0,2.0,3.0]
    node = ParameterElement('dvec', 'double',vec)
    input.attach(node)

    a = 3
    node = ParameterElement('ivalue','int',a) 

    # Branch Elements
    mesh = BranchElement('mesh')
    input.attach(mesh)
    cart_mesh = BranchElement('cartesian mesh')
    mesh.append(cart_mesh)

    base_num_blocks = [0,1]
    node = ParameterElement('base_num_blocks','int',base_num_blocks)
    cart_mesh.append(node)

    base_lower_left = [0.0,0.0]
    node = ParameterElement('base_lower_left','int',base_lower_left)
    cart_mesh.append(node)

    input.dumpXML()



