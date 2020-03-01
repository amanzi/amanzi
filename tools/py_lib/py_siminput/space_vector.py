#
# Copyright 2010 Los Alamos National Security, LLC
# Written by Robert B. Lowrie (CCS-2)
#
'''
Defines spatial vectors.
'''

from py_siminput.vector import Fixed_Vector

########################################################################

class Space_Vector(Fixed_Vector):
    '''
    General space vector.

    spec = initial vector.
    '''
    def __init__(self, spec):
        if len(spec) not in (1,2,3):
            raise TypeError('Initializer must be a list of length'\
                  ' 1,2, or 3')
        Fixed_Vector.__init__(self, 'd', spec)
        self._flag = 'vector_lite'
        
# Test:
if __name__ == '__main__':
    x = Space_Vector([0.2])
    print(x)

    x = Space_Vector([0.2, 33.4])
    print(x)

    x[1] = 2.2
    print(x)
