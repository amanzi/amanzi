#
# Copyright 2010 Los Alamos National Security, LLC
# Written by Robert B. Lowrie (CCS-2)
#
'''
Defines a vector class.
'''

from py_siminput.bool import Bool
from array            import array

class Vector:
    '''
    A vector class.  This class is similar to array.array, but
    also supports arrays of strings and bools.

    typecode : The type the Vector will store.  The types supported
               are as follows, with the default value given in brackets:
                 'd' for double     [0.0].
                 'i' for integer    [0]
                 'b' for bool.Bool  [Bool('false')]
                 's' for string     ['']
    spec     : Specifies the vector.  If spec a list of values, this list
               is used to initialize the vector.  The type of all the values in
               the list must correspond to one of the typecodes.  If spec is a
               non-negative integer, the vector length is set to spec and the
               element values are initialized to an appropriate default value.
    '''
    def __init__(self, typecode, spec):
        self._supported_types = ['d', 'i', 'b', 's']
        if typecode not in self._supported_types:
            raise ValueError('typecode must be one of ' + \
                  repr(self._supported_types))
        self.typecode = typecode
        
        # Create a list of values
        list = []
        if type(spec) == type(1):
            # if spec is an integer, create a list, of length spec,
            # using the default value to initialize.
            if spec < 0:
                raise ValueError('spec must be a list or a positive integer.')
            if spec > 0:
                defaults = {}
                defaults['d'] = 0.0;
                defaults['i'] = 0;
                defaults['b'] = Bool('false')
                defaults['s'] = ''
                list = spec * [defaults[typecode]]
        elif type(spec) == type([]):
            # spec is a list, so use it.
            list = spec
        else:
            raise ValueError('spec must be a list or an integer.')
        
        # Set the internal data
        if typecode in ['d', 'i']:
            # Use array for doubles and ints, because it handles various
            # numerical conversions nicely.
            self._use_array = 1
            self._data = array(typecode, list)
        else:
            # For bools and strings, keep track of the type ourselves
            # in attribute _type.
            self._use_array = 0
            if typecode == 'b':
                self._type = type(Bool())
            elif typecode == 's':
                self._type = type('a string')
                for i in list:
                    if type(i) != self._type:
                        raise TypeError(repr(i) + ' is illegal type.')
            self._data = list
        # This flag may be used by derived classes (see space_vector)
        self._flag = "std"
    def _check_type(self, value):
        '''
        Ensures that value is the correct type.
        '''
        # If using an array, the array type will do its own checking.
        if not self._use_array:
            if type(value) != self._type:
                raise TypeError(repr(value) + ' is illegal type.  ' + \
                      'Must be type ' + repr(self._type))
    def __len__(self):
        return len(self._data)
    def __getitem__(self, key):
        return self._data[key]
    def __setitem__(self, key, value):
        self._check_type(value)
        self._data[key] = value
    def __repr__(self):
        if self._use_array:
            return repr(self._data.tolist())
        else:
            return repr(self._data)
    def append(self, value):
        '''
        Appends value onto the end of the vector.
        '''
        self._check_type(value)
        self._data.append(value)


class Fixed_Vector(Vector):
    '''
    A fixed_vector is a vector of fixed length.  This is enforced by
    type_check.
    '''
    def __init__(self, typecode, spec):
        Vector.__init__(self, typecode, spec)

# Test code:
if __name__ == '__main__':
    x = Vector('d', [0.])
    x[0] = 1.0
    x[0] = 2
    x.append(2.7)
    print(repr(x))
    s = Vector('s', ['cat', 'dog'])
    s[1] = 'tiger'
    print(repr(s))
    b = Vector('b', [Bool('false'), Bool('true')])
    b.append(Bool('false'))
    print(repr(b))
    print(len(b))
    q = Vector('d', 3)
    print(repr(q))

