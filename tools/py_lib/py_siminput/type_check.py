#
# Copyright 2010 Los Alamos National Security, LLC
# Written by Robert B. Lowrie (CCS-2)
#
'''
Various tools for type-checking.
'''

import types
import py_siminput.vector
    
########################################################################

def create_enum(list, firstIndex = 0):
    '''
    Returns a dictionary of list:index.  This can be used for creating
    an enumerated-like type from a list of strings, to be passed as
    the dict option in the Type_Check._in_dict() function.
    '''
    i = firstIndex
    dict = {}
    for key in list:
        dict[key] = i
        i = i + 1
    return dict

########################################################################

def add_to_enum(enum, list):
    'Adds keys in list to enum.'
    i = max(enum.values()) + 1
    dict = enum.copy()
    for key in list:
        dict[key] = i
        i = i + 1
    return dict

########################################################################

def get_type(value):
    '''
    Extension of type builtin, except that:
       - for InstanceType, returns class type.
       - for Vector, adds the typecode
    '''
    tvalue = type(value)
    if isinstance(tvalue,object):
        tvalue = str(value.__class__)
        if isinstance(value, py_siminput.vector.Fixed_Vector):
            tvalue = tvalue + '(' + value.typecode + str(len(value)) + ')'
        elif isinstance(value, py_siminput.vector.Vector):
            tvalue = tvalue + '(' + value.typecode + ')'
    return tvalue

########################################################################

def compatible_types(v, ref):
    '''
    Returns true if v and ref are compatible numeric types.
    '''
    same = 1
    if v != ref:
        allInt = [int]
        allNum = allInt + [float]
        if isinstance(ref,float) and v in allNum:
            pass
        elif isinstance(ref,int) and v in allInt:
            pass
        else:
            same = 0
    return same

########################################################################

class _Constraint:
    '''
    Base constraint for all contraints.  This constraint will enforce
    that values be of a certain type.

    Attributes:

    const     = If set, value cannot be changed from def_value.
    def_value = The default value.
    type      = The type of def_value, as returned by get_type().
    '''
    def __init__(self, def_value, const):
        self.const     = const
        self.def_value = def_value
        self.type      = get_type(def_value)
        self.dump_type = self.type
    def check_type(self, value, name):
        '''
        Checks whether value is compatible with the current type.
        '''
        tv = get_type(value)
        if not compatible_types(tv, self.type):
            raise TypeError(name + ' is ' + str(self.type) + \
                  '; incompatible with ' + str(tv))
    def doc_constraints(self, obj, name):
        'Forms a documentation string of the constraints.'
        doc  = ''
        also = 'Must satisfy '
        if self.is_const():
            doc = doc + name + \
                  ' is a constant value'
            also = ', '
        return (doc, also)
    def print_value_error(self, value, obj, name):
        'Raises ValueError and error message.'
        doc = name + ' = ' + str(value) + ' invalid.\n'
        doc = doc + '       ' + self.doc_constraints(obj, name)[0]
        raise ValueError(doc)
    def check(self, value, obj, name):
        'Checks all constraints on value.'
        self.check_type(value, name)
    def is_const(self):
        return self.const
    def dump_value(self, v):
        'Returns the value to dump, corresponding to v.'
        return v

########################################################################

class _Constraint_Numeric(_Constraint):
    '''
    Constraint wrapper for numeric values.  Adds to _Constraint the
    ability to set minimum and maximum values.

    Attributes (in addition to those in _Constraint):

    min       = If set, minimum allowable value.
    max       = If set, maximum allowable value.
    min_attr  = A list of attributes which the value must be greater than.
    max_attr  = A list of attributes which the value must be less than.
    '''
    def __init__(self, def_value, const):
        _Constraint.__init__(self, def_value, const)
        self.min = None
        self.max = None
        self.min_attr = []
        self.max_attr = []
    def doc_constraints(self, obj, name):
        'Forms a documentation string of the constraints.'
        (doc, also)  = _Constraint.doc_constraints(self, obj, name)
        if self.min is not None:
            doc = doc + also + name + ' >= ' + str(self.min)
            also = ', '
        for v in self.min_attr:
            doc = doc + also + name + ' >= ' + v + ' = ' + str(getattr(obj, v))
            also = ', '
        if self.max is not None:
            doc = doc + also + name + ' <= ' + str(self.max)
            also = ', '
        for v in self.max_attr:
            doc = doc + also + name + ' <= ' + v + ' = ' + str(getattr(obj, v))
            also = ', '
        return (doc, also)
    def check_min(self, value, min, obj, name):
        if value < min:
            self.print_value_error(value, obj, name)
    def check_max(self, value, max, obj, name):
        if value > max:
            self.print_value_error(value, obj, name)
    def check(self, value, obj, name):
        'Checks all constraints on value.'
        self.check_type(value, name)
        if self.min is not None:
            self.check_min(value, self.min, obj, name)
        for v in self.min_attr:
            self.check_min(value, getattr(obj, v), obj, name)
        if self.max is not None:
            self.check_max(value, self.max, obj, name)
        for v in self.max_attr:
            self.check_max(value, getattr(obj, v), obj, name)
    def add_min(self, minValue):
        'Adds a minimum value constraint.'
        if isinstance(minValue,str):
            if minValue not in self.min_attr:
                self.min_attr.append(minValue)
        else:
            self.check_type(minValue, 'minValue')
            self.min = minValue
    def add_max(self, maxValue):
        'Adds a maximum value constraint.'
        if isinstance(maxValue,str):
            if maxValue not in self.max_attr:
                self.max_attr.append(maxValue)
        else:
            self.check_type(maxValue, 'maxValue')
            self.max = maxValue

########################################################################

class _Constraint_Class(_Constraint):
    '''
    Constraint that a value must have a base class.

    Attributes:

    base_class = Values must have this class as a base class.
    '''
    def __init__(self, def_value, const, base_class):
        _Constraint.__init__(self, def_value, const)
        if not isinstance(def_value, base_class):
            raise TypeError('def_value is not an instance of base_class')
        self.base_class = base_class
    def check_type(self, value, name):
        '''
        Checks whether value is compatible with the current type.
        '''
        tv = get_type(value)
        # ensure value is an instance of base_class
        if not isinstance(value, self.base_class):
            raise TypeError(name + ' is an instance of ' + \
                  str(self.base_class) + \
                  '; incompatible with type ' + str(tv))
    def doc_constraints(self, obj, name):
        'Forms a documentation string of the constraints.'
        (doc, also)  = _Constraint.doc_constraints(self, obj, name)
        doc = doc + also + name + \
              ' must be derived from a base class ' + str(self.base_class)
        also = ', '
        return (doc, also)
    def check(self, value, obj, name):
        'Checks all constraints on value.'
        self.check_type(value, name)

########################################################################

class _Constraint_List(_Constraint):
    '''
    Constraint that a value must be one of a list.

    Attributes:

    list = The list the value must be in.
    '''
    def __init__(self, def_value, const, list):
        _Constraint.__init__(self, def_value, const)
        for v in list:
            self.check_type(v, 'value in list')
        self.list = list
    def doc_constraints(self, obj, name):
        'Forms a documentation string of the constraints.'
        (doc, also)  = _Constraint.doc_constraints(self, obj, name)
        doc = doc + name + \
              ' must be one of the following: ' + str(self.list)
        also = ', '
        return (doc, also)
    def check(self, value, obj, name):
        'Checks all constraints on value.'
        self.check_type(value, name)
        if value not in self.list:
            self.print_value_error(value, obj, name)

########################################################################

class _Constraint_Type_List(_Constraint):
    '''
    Constraint that the type of a value must be one of a list.

    Attributes:

    list = The list of allowable types.

    Note that self.type refers to the type of the def_value, not of
    the current value.
    '''
    def __init__(self, def_value, const, list):
        _Constraint.__init__(self, def_value, const)
        self.list = []
        for v in list:
            self.append(get_type(v))
    def append(self, type):
        self.list.append(type)
    def doc_constraints(self, obj, name):
        'Forms a documentation string of the constraints.'
        (doc, also)  = _Constraint.doc_constraints(self, obj, name)
        doc = doc + name + \
              ' must be a type of one of the following: ' + str(self.list)
        also = ', '
        return (doc, also)
    def check(self, value, obj, name):
        'Checks all constraints on value.'
        if get_type(value) not in self.list:
            self.print_value_error(value, obj, name)

########################################################################

class _Constraint_Dict(_Constraint):
    '''
    Constraint that a value must be a key in a dictionary.

    Attributes:

    dict = Value must be in dict.keys().
    '''
    def __init__(self, def_value, const, dict):
        _Constraint.__init__(self, def_value, const)
        for key in dict.keys():
            self.check_type(key, 'key in dict')
        self.dict = dict
        self.dump_type = get_type(dict[def_value])
    def doc_constraints(self, obj, name):
        'Forms a documentation string of the constraints.'
        (doc, also)  = _Constraint.doc_constraints(self, obj, name)
        doc = doc + name + \
              ' must be one of the following: ' + str(self.dict.keys())
        also = ', '
        return (doc, also)
    def check(self, value, obj, name):
        'Checks all constraints on value.'
        self.check_type(value, name)
        if value not in self.dict.keys():
            self.print_value_error(value, obj, name)
    def dump_value(self, v):
        '''
        Returns the value to dump, corresponding to v.  In this case,
        it is the corresponding token value.
        '''
        return self.dict[v]

########################################################################

class Type_Check:
    '''
    A class that typechecks its attributes.  This class controls
    assignment of its attributes, such as

        self.name = value

    The above operation is allowed only if:

        1) name has been declared.
        2) value is a compatible type with that defined when name was
           declared.
        3) value satisfies all constraints (min, max, etc.) set for
           name.

    name is declared using the member functions _declare, _numeric,
    _derived, _in_dict, or _in_list.  We prefix all functions with an
    underscore in order to avoid attribute name conflicts.
        
    Internal attribute:

    _constraint = Dictionary of name string to Constraint object.
    '''
    def __init__(self): 
        # dictionary of declared names to constraint objects
        self._set('_constraint', {})
    def _set(self, name, value):
        '''
        Does self.name = value with no checking.  Useful for internal
        attributes.
        '''
        self.__dict__[name] = value
    def _check_if_declared(self, name):
        '''
        Checks if name has been declared.  Throws AttributeError if not.
        '''
        if name not in self._constraint.keys():
            raise AttributeError(name + ' has not been declared')
    def _check_if_reserved(self, name):
        '''
        Checks if name is a reserved attribute.
        '''
        if hasattr(self, name) and name not in self._constraint.keys():
            raise AttributeError(name + \
                  ' is a reserved attribute name for ' + str(self.__class__))
    def _declare(self, name, value, const=False):
        '''
        Declare an attribute name and set its initial value.
        Later uses of

            self.name = new_value

        will enforce that new_value is the same type as value.

        name      = attribute string name.
        value     = initial value.
        const     = If true, value cannot be changed from initial value.
        '''
        self._check_if_reserved(name)
        self._constraint[name] = _Constraint(value, const)
        self._set(name, value)
    def _numeric(self, name, value, min=None, max=None, const=False):
        '''
        Declare a numeric attribute name and set its initial value.
        Later uses of

            self.name = new_value

        will enforce that new_value is a compatible numeric type as
        value (see compatible_types()) and that new_value satisfies
        min and max constraints.

        name      = attribute string name.
        value     = initial value.
        min       = If set, minimum value, or list of mins.
        max       = If set, maximum value, or list of maxs.
        const     = If True, value cannot be changed from initial value.
        '''
        self._check_if_reserved(name)
        c = _Constraint_Numeric(value, const)
        self._set(name, value)
        # Add min constraints
        if min is not None:
            if not isinstance(min,list):
                min = [min]
            for z in min:
                if isinstance(z,str):
                    self._check_if_declared(z)
                    self._constraint[z].add_max(name)
                c.add_min(z)
                c.check(getattr(self, name), self, name)
        # Add max constraints
        if max is not None:
            if not isinstance(max,list):
                max = [max]
            for z in max:
                if isinstance(z,str):
                    self._check_if_declared(z)
                    self._constraint[z].add_max(name)
                c.add_max(z)
                c.check(getattr(self, name), self, name)
        self._constraint[name] = c
    def _in_dict(self, name, value, dict, const=False):
        '''
        Declare an attribute name and set its initial value.
        Later uses of

            self.name = new_value

        will enforce that new_value is a key in dict.  Certain
        functions (such as Dump functions) dump the token of value,
        instead of value itself.

        name      = attribute string name.
        value     = initial value.
        dict      = value must be in keys of this dictionary.
        const     = If True, value cannot be changed from initial value.
        '''
        self._check_if_reserved(name)
        self._constraint[name] = _Constraint_Dict(value, const, dict)
        self._set(name, value)
    def _in_list(self, name, value, list, const=False):
        '''
        Declare an attribute name and set its initial value,
        constrained to a list of values.

        name      = attribute string name.
        value     = initial value.
        list      = value must be in this list.
        const     = If True, value cannot be changed from initial value.
        '''
        self._check_if_reserved(name)
        self._constraint[name] = _Constraint_List(value, const, list)
        self._set(name, value)
    def _type_list(self, name, value, list, const=False):
        '''
        Declare an attribute name and set its initial value,
        with its type constrained to a list of types.

        name      = attribute string name.
        value     = initial value.
        list      = the type of value must be in this list.
        const     = If True, value cannot be changed from initial value.
        '''
        self._check_if_reserved(name)
        self._constraint[name] = _Constraint_Type_List(value, const, list)
        self._set(name, value)
    def _derived(self, name, value, base_class, const=False):
        '''
        Declare an attribute name and set its initial value,
        constrained to have a base class of a given type.

        name       = attribute string name.
        value      = initial value.
        base_class = value must be or have a base class of this type.
        const      = If True, value cannot be changed from initial value.
        '''
        self._check_if_reserved(name)
        self._constraint[name] = _Constraint_Class(value, const, base_class)
        self._set(name, value)
    def _set_to_default(self, name = None):
        'Sets name (or all declared names if None) back to its default value.'
        if name is None:
            for key in self._constraint.keys():
                self._set_to_default(key)
        else:
            self._check_if_declared(name)
            v = getattr(self, name)
            if isinstance(v, Type_Check):
                v._set_to_default()
            else:
                self._set(name, self._constraint[name].def_value)
    def __setattr__(self, name, value):
        'Does self.name = value with some checking.'
        self._check_if_declared(name)
        if self._constraint[name].is_const():
            raise TypeError(name + ' is const.')
        self._constraint[name].check(value, self, name)
        self._set(name, value)
    def _dump_value(self, name):
        '''
        Returns the value to dump corresponding to the attribute name.
        Typically, this is self.name, but for Dict constraints, it is
        the corresponding token value.
        '''
        self._check_if_declared(name)
        return self._constraint[name].dump_value(getattr(self, name))
