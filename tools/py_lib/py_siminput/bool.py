#
#  $Id: bool.py 724 2009-07-22 19:17:03Z wohlbier $
#
'''
Defines a boolean class.

Written by Rob Lowrie, May 2004
'''

class Bool:
    '''
    A boolean class.
    '''
    def __init__(self, b = None):
        self.false  = ['false', 'f', 0, False]
        self.true   = ['true', 't', 1, True]
        self.state = 0 # default to false
        if b is not None:
            self._set(b)
    def __call__(self, b = None):
        if b is not None:
            self._set(b)
        return self.state
    def __repr__(self):
        if self.state:
            return 'true'
        else:
            return 'false'
    def _set(self, b):
        if b in self.true:
            self.state = 1
        elif b in self.false:
            self.state = 0
        else:
            raise TypeError(b + " unrecogonized value for Bool.")

# Test code:
if __name__ == '__main__':
    x= Bool()
    if ( x() ): print('failed.')
    if ( x(0) ): print('failed.')
    if ( x('f') ): print('failed.')
    if ( x('false') ): print('failed.')
    if ( x(False) ): print('failed.')
    if ( not x(1) ): print('failed.')
    if ( not x('t') ): print('failed.')
    if ( not x('true') ): print('failed.')
    if ( not x(True) ): print('failed.')
    print(repr(x), str(x))
    x = Bool('f')
    print(repr(x), str(x))

