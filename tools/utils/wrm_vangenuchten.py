import numpy as np

class VanGenuchten(object):
    def __init__( self, m=0., alpha=0., sr=0.0 ):
        self._m = m
        self._alpha = alpha
        self._sr = sr
        self._n = 1.0/(1.0-m)

    def capillaryPressure( self, s ):
        se = (s - self._sr) / (1.0 - self._sr);
        return (pow(pow(se, -1.0/self._m) - 1.0, 1/self._n)) / self._alpha

    def d_saturation( self, pc ):
        if pc > 0.0:
            return -self._m*self._n * pow(1.0 + pow(self._alpha*pc, self._n), -self._m-1.0) * pow(self._alpha*pc, self._n-1) * self._alpha * (1.0 - self._sr)
        else:
            return 0.0

    def saturation( self, pc ):
        if (pc > 0.0):
            return pow(1.0 + pow(self._alpha*pc, self._n), -self._m) * (1.0 - self._sr) + self._sr
        else:
            return 1.0

    def k_relative( self, pc ):
        se = pow(1.0 + pow(self._alpha*pc,self._n),-self._m)
        return np.sqrt(se) * pow( 1.0 - pow( 1.0 - pow(se,1.0/self._m),self._m), 2);

