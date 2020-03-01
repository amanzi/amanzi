#
# Copyright 2010 Los Alamos National Security, LLC
# Written by Robert B. Lowrie (CCS-2)
#
'''
Interface for user code input.
'''

import sys, string, types
import xml.sax.saxutils

from py_siminput.type_check   import Type_Check, create_enum, get_type
from py_siminput.bool         import Bool
from py_siminput.vector       import Vector
from py_siminput.space_vector import Space_Vector

########################################################################

def format_paragraph(s, indent='', width=75):
    '''
    Formats a string into a paragraph by inserting newline characters.

    s      = The input string.
    indent = A string added to the beginning of each line.
    width  = The maximum width of a line (including indent).

    Returns: A string, with newline characters and indent strings inserted.
    '''
    p = indent
    nIndent = len(indent)
    currentWidth = nIndent
    for word in string.split(s):
        if currentWidth > width:
            currentWidth = nIndent
            p= p[0:len(p)-1]  + '\n' + (nIndent * ' ')
        add = word + ' '
        p = p + add
        currentWidth = currentWidth + len(add)
    p = p + '\n'
    return p

########################################################################

class Dump_XML:
    '''
    Dumps an XML file of an Input_Interface.

    Attributes:

      indent_increment = The increment to indent each level of nested
                         entities.
      type_map         = dictionary that maps types to strings.
    '''
    def __init__(self, x = None, filename = None):
        '''
        Constructor.

        x = If not None, calls dump(x, filename).
        '''
        # the amount to incremented nested xml entities
        self.indent_increment = '   '
        # define the mapping of types
        self.type_map = {}
        self.type_map[float]  = 'real'
        self.type_map[int]    = 'int'
        self.type_map[int]   = 'int'
        self.type_map[bytes] = 'string'
        self.type_map[get_type(Bool())]  = 'bool'
        if x is not None:
            self.dump(x, filename)
    def dump(self, x, filename = None):
        '''
        Outputs an Input_Interface declared attributes in an XML
        format.

        x        = A InterInterface object.
        filename = file name to dump output to.  If None, stdout.
        '''
        if filename is None:
            file = sys.stdout
        else:
            file = open(filename, 'w')
        file.write('<?xml version="1.0" encoding="iso-8859-1"?>\n')
        file.write('<!DOCTYPE branch [\n')
        file.write('   <!ELEMENT branch (value*)>\n')
        file.write('   <!ATTLIST branch label CDATA #REQUIRED>\n')
        file.write('   <!ELEMENT value (bool|branch|int|real|string|vector)+>\n')
        file.write('   <!ATTLIST value label CDATA #REQUIRED>\n')
        file.write('   <!ELEMENT bool (#PCDATA)>\n')
        file.write('   <!ELEMENT int (#PCDATA)>\n')
        file.write('   <!ELEMENT real (#PCDATA)>\n')
        file.write('   <!ELEMENT string (#PCDATA)>\n')
        file.write('   <!ELEMENT vector (entry*)>\n')
        file.write('   <!ATTLIST vector type (bool|int|real|string) #REQUIRED\n')
        file.write('                    length CDATA #REQUIRED\n')
        file.write('                    flag CDATA #REQUIRED\n')
        file.write('                    default CDATA #IMPLIED>\n')
        file.write('   <!ELEMENT entry (#PCDATA)>\n')
        file.write('   <!ATTLIST entry index CDATA #REQUIRED>\n')
        file.write(']>\n')
        self._dump_recursive(x, '', file)
        if file != sys.stdout:
            file.close()
    def _dump_recursive(self, x, indent, file):
        '''
        Used by dumpXML to recursively dump input_interface entries.
        '''
        file.write(indent + '<branch label="' + x._label + '">\n')
        ip  = indent + self.indent_increment
        ip2 = ip + self.indent_increment
        for name in x._print_order:
            file.write(ip + '<value label="' + name + '">\n')
            v = x._dump_value(name)
            if isinstance(v, Input_Interface):
                # Allow nested Input_Interfaces
                self._dump_recursive(v, ip2, file)
            elif isinstance(v, Vector):
                # We have a vector
                self._write_vector(v, ip2, file)
            else:
                vtype = x._constraint[name].dump_type
                self._write_std(v, vtype, ip2, file)
            file.write(ip + '</value>\n')
        file.write(indent + '</branch>\n')
    def _write_vector(self, v, indent, file):
        '''
        Writes a vector in XML.

        v      = an array.
        indent = current indentation level.
        file   = an open file stream.
        '''
        vector_map = {}
        vector_map['d']  = 'real'
        vector_map['i']  = 'int'
        vector_map['s']  = 'string'
        vector_map['b']  = 'bool'
        file.write(indent + '<vector')
        file.write(' type="' + vector_map[v.typecode] + '"')
        file.write(' length="' + repr(len(v)) + '"')
        file.write(' flag="' + v._flag + '"')
        file.write('>\n')
        ip = indent + self.indent_increment
        for i in range(len(v)):
            file.write(ip + \
                       '<entry index="%s">%s</entry>\n' \
                       % (i, str(v[i])))
        file.write(indent + '</vector>\n')
    def _write_std(self, v, vtype, indent, file):
        '''
        Writes standard types in XML.  Standard types are those
        defined in the type_map attribute.

        v      = A value
        vtype  = The type of v.  We pass this in because we want to match
                 the type of the original declaration (so for example,
                 if v is type int, but declared type real).
        indent = Current indentation level.
        file   = An open file stream.
        '''
        if vtype not in list(self.type_map.keys()):
            raise TypeError(repr(v) + ' has unrecogonized type: ' + \
                  repr(vtype) + '\nMust be one of ' + \
                  repr(list(self.type_map.keys())))
        mtype = self.type_map[vtype]
        file.write(indent + '<' + mtype + '>')
        file.write(xml.sax.saxutils.escape(str(v)))
        file.write('</' + mtype + '>\n')

########################################################################

class Dump_Output:
    '''
    Dumps an Input_Interface as human-readable output.  The
    human-readable format is basically a flat file, with one value per
    line.
    '''
    def __init__(self, x = None, filename = None):
        '''
        Outputs an Input_Interface declared attributes in a human-readable
        format.

        x = If not None, dump(x, filename) is called.
        '''
        if x is not None:
            self.dump(x, filename)
    def dump(self, x, filename = None):
        '''
        Outputs an Input_Interface declared attributes in a human-readable
        format.

        x        = A InterInterface object.
        filename = file name to dump output to.  If None, stdout.
        '''
        if filename is None:
            file = sys.stdout
        else:
            file = open(filename, 'w')
        self._dump_recursive(x, file)
        if file != sys.stdout:
            file.close()
    def _dump_recursive(self, x, file, parents = None):
        '''
        Recursive dump function.

        x       = A InterInterface object.
        file    = An open file stream.
        parents = A concatenation of parent, grandparent, ... variable
                  names, for example great_grandparent.grandparent.parent
        '''
        for name in x._print_order:
            if parents is None:
                fullName = name
            else:
                fullName = parents + '.' + name
            v = x._dump_value(name)
#            if a is not None:
#                # We have a dictionary
#                fullName = '(' + `v` + ') ' + fullName
            if isinstance(v, Input_Interface):
                # Allow nested Input_Interfaces
                self._dump_recursive(v, file, fullName)
            elif isinstance(v, Vector):
                # We have a vector
                file.write(repr(len(v)) + ' # ' + fullName + '\n')
                for i in range(len(v)):
                    file.write(str(v[i]) + '\n')
            else:
                file.write(repr(v) + ' # ' + fullName + '\n')

########################################################################

class Dump_Python:
    '''
    Dumps an Input_Interface as python input.

    The output from this function will most likely need to be edited,
    in order to use it as actual input.  In particular:

       - The needed modules must be imported.
       - Objects are created assuming that they have a constructor
       (__init__) that takes no arguments.  This may not be the case.

    Attributes:

      doc_comments = If true, output documentation as comments before
                     each variable.
    '''
    def __init__(self, doc_comments = 0):
        self.doc_comments = doc_comments
    def dump(self, x, var_name, filename = None):
        '''
        Outputs Input_Interface declared attributes as python input.

        x        = A InterInterface object.
        var_name = The variable name to be used for x.
        filename = file name to dump output to.  If None, stdout.
        '''
        if filename is None:
            file = sys.stdout
        else:
            file = open(filename, 'w')
        self._dump_recursive(x, var_name, file)
        if file != sys.stdout:
            file.close()
    def _dump_recursive(self, x, var_name, file):
        for name in x._print_order:
            if self.doc_comments:
                file.write('\n')
                if x._doc[name] != 'NO':
                    file.write(format_paragraph(x._doc[name], '# '))
            fullName = var_name + '.' + name
            v = getattr(x, name)
            if isinstance(v, Input_Interface):
                # Allow nested Input_Interfaces.
                #
                # This write assumes that the object's __init__ function
                # takes no arugments.  It's difficult to do much better
                # than this.
                file.write(fullName + ' = ' + str(v.__class__) + \
                           '()\n')
                self._dump_recursive(v, fullName, file)
            elif isinstance(v, Space_Vector):
                # We have a space vector
                file.write(fullName + ' = Space_Vector([')
                for i in range(len(v)):
                    if i > 0: file.write(',')
                    file.write(repr(v[i]))
                file.write(']\n')
            elif isinstance(v, Vector):
                # We have a vector
                file.write(fullName + ' = Vector(' + \
                           repr(v.typecode) + ', ' + repr(len(v)) + ')\n')
                for i in range(len(v)):
                    file.write(fullName + '[' + repr(i) + \
                               '] = ' + repr(v[i]) + '\n')
            elif isinstance(v, Bool):
                # We have a bool
                file.write(fullName + ' = Bool("' + repr(v) + '")\n')
            else:
                file.write(fullName + ' = ' + repr(v) + '\n')

########################################################################

class Dump_Doc:
    '''
    Dumps the documentation of an Input_Interface.
    '''
    def __init__(self, x = None, filename = None):
        if x is not None:
            self.dump(x, filename)
    def dump(self, x, filename = None):
        '''
        Outputs the documenation for an Input_Interface.

        x        = A Input_Interface object.
        filename = file name to dump output to.  If None, stdout.
        '''
        if filename is None:
            file = sys.stdout
        else:
            file = open(filename, 'w')
        self._dump_recursive(x, file)
        if file != sys.stdout:
            file.close()
    def _dump_recursive(self, x, file, parents = None, basePad = ''):
        '''
        Recursive dump function.

        x       = A InterInterface object.
        parents = A concatenation of parent, grandparent, ... variable
                  names, for example great_grandparent.grandparent.parent
        basePad = padding string to start each line.
        '''
        for name in x._print_order:
            if x._doc[name] != 'NO':
                pad = basePad + '.... '
                file.write('\n' + basePad + name)
                if parents is None:
                    fullName = name
                else:
                    fullName = parents + '.' + name
                    file.write(' (' + parents + ')')
                file.write(' :\n')
                v = getattr(x, name)
                c = x._constraint[name]
                file.write(format_paragraph(x._doc[name], pad))
                # Use get_type(v) instead of c.type
                file.write(pad + repr(get_type(v)) + '\n')
                (doc, also) = c.doc_constraints(x, name)
                if doc:
                    file.write(format_paragraph(doc, pad))
                if isinstance(v, Input_Interface):
                    # Allow nested Input_Interfaces
                    self._dump_recursive(v, file, fullName, basePad + '     ')
                else:
                    file.write(pad + 'Present Value =  ' + repr(v) + '\n')

########################################################################

class Input_Interface(Type_Check):
    '''
    Input interface for simulation code input.  The attributes of this
    class are type checked by inheriting the Type_Check class.  This
    class adds the following internal attributes to Type_Check:

    _print_order = A list of declared attributes in the order they are to be
                   printed by the Dump classes.
    _label       = A string label (or name) of this object.
    _doc         = A dictionary of attribute name to documentation string.

    An Input_Interface is an n-tree, where at any level, attributes may
    include numeric data types, strings, dictionaries, or Input_Interfaces.
    '''
    def __init__(self): 
        Type_Check.__init__(self)
        self._set('_print_order', [])
        self._set('_label', 'none')
        self._set('_doc', {})
    def _set_doc(self, name, after_name, doc):
        '''
        Sets the documentation of name.
        '''
        # Only append to print_order if name is a new attribute
        if not hasattr(self, name):
            if after_name is None:
                self._print_order.append(name)
            else:
                i = self._print_order.index(after_name)
                self._print_order.insert(i+1, name)
        self._doc[name] = doc
    def _declare(self, name, value, const=False, after_name=None,
                 doc='UNDOCUMENTED'):
        '''
        See Type_Check._declare for more information.  This
        function adds the following parameters:

        after_name = If set, insert name in _print_order after this
                     name.
        doc        = Documentation string.  If NO, Dump_Doc does not
                     output documentation.
        '''
        self._set_doc(name, after_name, doc)
        Type_Check._declare(self, name, value, const)
    def _numeric(self, name, value, min=None, max=None, const=False,
                after_name=None, doc='UNDOCUMENTED'):
        '''
        See _declare above and Type_Check._numeric for documentation.
        '''
        self._set_doc(name, after_name, doc)
        Type_Check._numeric(self, name, value, min, max, const)
    def _in_dict(self, name, value, dict, const=False,
               after_name=None, doc='UNDOCUMENTED'):
        '''
        See _declare above and Type_Check._in_dict for documentation.
        '''
        self._set_doc(name, after_name, doc)
        Type_Check._in_dict(self, name, value, dict, const)
    def _in_list(self, name, value, list, const=False, after_name=None,
                doc='UNDOCUMENTED'):
        '''
        See _declare above and Type_Check._in_list for documentation.
        '''
        self._set_doc(name, after_name, doc)
        Type_Check._in_list(self, name, value, list, const)
    def _type_list(self, name, value, list, const=False, after_name=None,
                   doc='UNDOCUMENTED'):
        '''
        See _declare above and Type_Check._type_list for documentation.
        '''
        self._set_doc(name, after_name, doc)
        Type_Check._type_list(self, name, value, list, const)
    def _derived(self, name, value, base_class, const=False, after_name=None,
                doc='UNDOCUMENTED'):
        '''
        See _declare above and Type_Check._derived for documentation.
        '''
        self._set_doc(name, after_name, doc)
        Type_Check._derived(self, name, value, base_class, const)
