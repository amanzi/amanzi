################################################################################

import os, sys
import types

from xml.etree.ElementTree import Comment, ProcessingInstruction, QName
from xml.etree.ElementTree import _encode, _escape_cdata, _escape_attrib
from xml.etree.ElementTree import ElementTree, Element, _ElementInterface

################################################################################
'''
Global defines
'''
_ParameterListTag = 'ParameterList'
_ParameterTag     = 'Parameter'

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

def isparameter(object):
    return isinstance(object,_ParameterInterface)

def islist(object):
    return isinstance(object,ParameterList)

def islistroot(object):
    return isinstance(object,_ParameterListRootInterface)

def get_list_type(a):
    ret_type = type(a[0])
    for item in a:
        if type(item) != ret_type:
            ret_type = str
            break
    return ret_type

def get_str_list_type(a):

    tvalue = get_list_type(a)
    str_type = None
    
    if isinstance(tvalue,float):
        str_type = 'Array double'
    elif isinstance(tvalue,int):
        str_type = 'Array int'
    elif isinstance(tvalue,long):
        str_type = 'Array long'
    elif isinstance(tvalue,str):
        str_type = 'Array string'
    elif isinstance(tvalue,bool):
        str_type = 'Array bool'
    else:
        raise ValueError(str(value) + 'unknown value type' + str(tvalue))
    return str_type

def get_str_type(value):
    tvalue = type(value)
    str_type = None
    if isinstance(tvalue,float):
        str_type = 'double'
    elif isinstance(tvalue,int):
        str_type = 'int'
    elif isinstance(tvalue,long):
        str_type = 'long'
    elif isinstance(tvalue,str):
        str_type = 'string'
    elif isinstance(tvalue,bool):
        str_type = 'bool'
    elif isinstance(tvalue,object):
        str_type = str(value.__class__)
    elif isinstance(tvalue,list):
        str_type = get_str_list_type(value)
    else:
        raise ValueError(str(value) + 'unknown value type' + \
               str(tvalue))
    return str_type

def get_py_type(str_type):

    result = type(None)

    bool_type_strings = ['bool', 'boolean']
    int_type_strings  = ['int', 'integer']
    float_type_strings = ['real', 'double', 'float']
    string_type_strings = ['str', 'string', 'char', 'character']

    if str_type in bool_type_strings:
        result = bool
    elif str_type in int_type_strings:
        result = int
    elif str_type in float_type_strings:
        result = float
    elif str_type in string_type_strings:
        result = str
    else:
        raise ValueError(str_type + 'is an unknown string type')
     
    return result 


def convert_str_to_type(str,str_type):
 
    if not isinstance(str_type,str):
        raise TypeError(type(str_type).__name__ + ' is an invalid type')

    result = None
    py_type = get_py_type(str_type)
    if py_type is not None:
        if isinstance(py_type,int):
            result = int(str)
        elif isinstance(py_type,float):
            result = float(str)
        elif isinstance(py_type,bool):
            result = bool(str)
        elif isinstance(py_type,str):
            result = str
        else:
            raise ValueError(py_type + ' no conversion possible with this type')
    else:
        raise ValueError(str_type + ' can not convert this to a Python type')

    return result

################################################################################
#class TrilinosParser(XMLParser):
#
#   def __init__(self):
#
#       raise NotImplementedError, 'Trilinos Parser class is not implemented'
#
################################################################################
class _ParameterInterface(_ElementInterface):

    def __init__(self, name=None, value=None):
        
        if name is None:
            raise ValueError('Parameter constructor requires a name')

        if value is None:
            raise ValueError('Parameter constructor requires a value')
      
        #print "name=",name 
        #print "value=",value 
	#print "type=",get_str_type(value)

	attrib= {'name':name,'value':str(value),'type':get_str_type(value)}
        _ElementInterface.__init__(self,_ParameterTag,attrib)

    def get_name(self):
        return self.get('name',default=None)

    def set_name(self,name):
        return self.set('name',name)

    def get_type(self):
        str_type = self.get('type',default=None)
        return get_py_type(str_type)
    
    def get_value(self):
        value = self.get('value',default=None)
        str_type = self.get('type',default=None)
        return convert_str_to_type(value,str_type)

    def set_value(self,value):
        import re
        if isinstance(value,list):
            str1 = str(value)
            str2 = re.sub('\[','{',str1)
            str_value = re.sub('\]','}',str2)
        else:
            str_value = str(value)

        str_type = get_str_type(value)
        self.set('type', str_type)
        self.set('value', str_value)
        return str_value



################################################################################
class _ParameterListRootInterface(_ElementInterface):

    def __init__(self,name=None):

        if name == None:
            raise ValueError('Parameter constructor requires a name')

        _ElementInterface.__init__(self,_ParameterListTag,{})
        self.set('name',name)

    def get_name(self):
        return self.get('name')

################################################################################

def Parameter(name,value):
    return _ParameterInterface(name,value)

def ParameterListRoot(name):
    return _ParameterListRootInterface(name)

################################################################################
class ParameterList(ElementTree):

    def __init__(self,name=None, file=None):

        # Throw an error if name or file is not defined
        if name == None and file == None:
            raise ValueError('creating ParameterList instance requires a name' \
                              + ' or file')

        if name != None:
            root = ParameterListRoot(name)
	    ElementTree.__init__(self,element=root,file=None)
	else:
            ElementTree.__init__(self,element=None,file=file)

    def name(self):
        root = self.getroot()
        return root.get('name')

    def attach(self,object):
        root = self.getroot()
        if isparameter(object) or islistroot(object):
            root.append(object)
        elif islist(object):
            oroot = object.getroot()
            root.append(oroot) 
        else:
            raise TypeError('Can not attach type =' + get_str_type(object) + \
                  ' is not a valid parameter or list element')

    def add_parameter(self,name,data):
        node = Parameter(name,data)
        self.attach(node)
        return node

    def add_sublist(self,object):
        if isinstance(object,str):
            list = ParameterList(object)
            root = list.getroot()
            self.attach(root)
            return list
        elif islist(object):
            root = object.getroot()
            self.attach(root)
            return object
        else:
            raise TypeError('type=' + type(object).__name__ \
                              + ' is an invalid type')

    def find_sublist(self,target):
        sublist = None
        result = None
        search_list = []
        try:
            search_list = self.iterfind(_ParameterListTag)
        except AttributeError:
            search_list = self.getiterator(tag=_ParameterListTag)

        for node in search_list:
            node_name = node.get('name')
            if node_name == target:
                result = node
                break
                 
        if result == None:
            print('Could not find sublist with name=' + target)
        else:
            list_name = result.get('name')
            sublist = ParameterList(list_name)
            sublist._setroot(result)
         
        return sublist

    def find_parameter(self,target=None,quiet=False):
        result = None

        if target == None:
            raise ValueError(' requires a parameter name')

	'''
        search_list = []    
        try:
            search_list = self.iterfind(_ParameterTag)
        except AttributeError:
            search_list = self.getiterator(_ParameterTag)
	'''
        
	search_list = self.findall(_ParameterTag)
        index=0
        for node in search_list:
            node_name = node.get('name')
            if node_name == target:
	        tri_type=node.get('type')
		py_value=convert_str_to_type(node.get('value'),tri_type) 
	        result = Parameter(node.get('name'),py_value)
		root = self.getroot()
		try:
		  root.remove(node)
		except ValueError:
		  print(root.tag)
		  print(root.get('name'))
		  for item in root.getchildren():
		    print("CHILD tag", item.tag, " name=", item.get('name'))
		  raise
		root.insert(index,result)
                break
	    index=index+1  

        if result == None and quiet != False:
            print('Could not find parameter ' + target)

        return result    

    def set_parameter(self,name,value):
        node = self.find_parameter(name,quiet=True)
        if node != None:
            node.set_value(value)
        else:
            node = Parameter(name,value)
            self.attach(node)

        return node

    def add_verbose(self,level=None):
	verbose = self.add_sublist("VerbosityObject")
	verbose.add_parameter("Verbosity Level", str(level))
	return verbose

    def dumpXML(self,file=sys.stdout,encoding='utf-8',xml_translate=None,*args):

        sortflag='name,type,value'
        if xml_translate != None :
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
                        except:
                            sys.stderr.write("sortcmp not a valid comparator, sorting alphabetically instead")
                            items.sort()
                else:
                    try:
                        items.sort(cmp=sortcmp)
                    except:
                        sys.stderr.write("sortcmp not a valid comparator, sorting alphabetically instead")
                        items.sort()
                ###

                        
                for k, v in items:
                    try:
                        if isinstance(k, QName) or k[:1] == "{":
                            k, xmlns = fixtag(k, namespaces)
                            if xmlns: xmlns_items.append(xmlns)
                    except TypeError:
                        _raise_serialization_error(k)
                    try:
                        if isinstance(v, QName):
                            v, xmlns = fixtag(v, namespaces)
                            if xmlns: xmlns_items.append(xmlns)
                    except TypeError:
                        _raise_serialization_error(v)
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
def InputList(file=None):

    if file == None:
        return ParameterList('Main')
    else:
        return ParameterList(file=file)

################################################################################
        # If run as a script, do some testing
if __name__ == '__main__':

   # Set up some input parameters and dump to STDOUT
   input = InputList()

   # Sample mesh
   mesh = ParameterList('Mesh')
   mesh.add_parameter('Mesh Class','MOAB')
   moab_list = mesh.add_sublist('MOAB Mesh Parameters')
   moab_list.add_parameter('Exodus file name','fbasin_unstr_400_V02.exo')

   # MPC
   mpc = ParameterList('MPC')
   mpc.add_parameter('Start Time', 0.0)
   mpc.add_parameter('End Time', 315360000.0)
   mpc.add_parameter('End Cycle', -1)
   mpc.add_parameter('disable Flow PK','no')
   mpc.add_parameter('disable Transport PK','no')
   mpc.add_parameter('disable Chemistry PK','no')
   mpc.add_parameter('Viz dump cycle frequency', -1)
   mpc.add_parameter('Viz dump time frequency',2592000.0 )
   cgns = mpc.add_sublist('CGNS')
   cgns.add_parameter('File name','fbasin-5-component.cgns')

   input.add_sublist(mesh)
   input.add_sublist(mpc)

   input.dumpXML()

   mpc.set_parameter('End Cycle', 100)
   mpc.dumpXML()


   # Example of an array parameter
   array_list = ParameterList("Array List")
   a = [0.0, 0.1, 0.2]
   array_list.add_parameter("Double Array",a)
   a = [0, 1, 2]
   array_list.add_parameter("Int Array",a)
   array_list.dumpXML()


   # Read Fbasin input file
   #fbasin = InputList(file='fbasin-5-components-025.xml')
   #print type(fbasin).__name__
   #print fbasin
   #flow = fbasin.find_sublist('Flow')
