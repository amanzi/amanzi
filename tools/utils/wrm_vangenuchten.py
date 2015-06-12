import numpy as np

class VanGenuchten(object):
    def __init__( self, m=0., alpha=0., sr=0.0, smoothing_interval=0. ):
        self._m = m
        self._alpha = alpha
        self._sr = sr
        self._n = 1.0/(1.0-m)
        self._pc0 = smoothing_interval

        if self._pc0 > 0.:
            k0 = self.k_relative(self._pc0) - 1.
            k0p = self.d_k_relative(self._pc0)
            self._a = (3 * k0 - k0p * self._pc0) / (self._pc0**2)
            self._b = (k0p * self._pc0 - 2 * k0) / (self._pc0**3)

    def capillaryPressure( self, s ):
        se = (s - self._sr) / (1.0 - self._sr);
        if (se < 1.e-8):
            return pow(se, -1.0/(self._m * self._n)) / self._alpha
        else:
            return (pow(pow(se, -1.0/self._m) - 1.0, 1/self._n)) / self._alpha

    def d_capillaryPressure( self, s ):
        se = (s - self._sr) / (1.0 - self._sr);
        if se < 1.e-8:
            return -1 / (self._m * self._n * self._alpha) \
              * pow(se, -1/(self._m * self._n) - 1.) / (1. - self._sr)
        else:
            return -1 / (self._m * self._n * self._alpha) \
              * pow(pow(se,-1./self._m) - 1., 1./self._n - 1.) \
              * pow(se, -1./self._m - 1.) / (1. - self._sr)

    def saturation( self, pc ):
        if pc <= 0.0:
            return 1.0
        else:
            se = pow(1.0 + pow(self._alpha*pc, self._n), -self._m)
            return  se * (1.0 - self._sr) + self._sr

    def d_saturation( self, pc ):
        if pc <= 0.:
            return 0.
        else:
            dse_dpc = -self._m*self._n * pow(1.0 + pow(self._alpha*pc, self._n), -self._m-1.0) * pow(self._alpha*pc, self._n-1) * self._alpha
            return  (1.0 - self._sr) * dse_dpc

    def k_relative( self, pc ):
        if pc >= self._pc0:
            se = pow(1.0 + pow(self._alpha*pc,self._n),-self._m)
            return np.sqrt(se) * pow( 1.0 - pow( 1.0 - pow(se,1.0/self._m),self._m), 2);
        elif pc <= 0.:
            return 1.
        else:
            return 1. + self._a * (pc**2) + self._b * (pc**3)

    def d_k_relative( self, pc ):
        if pc >= self._pc0:
            se = pow(1.0 + pow(self._alpha*pc,self._n),-self._m)
            dsdp = self.d_saturation(pc)

            x = pow(se,1./self._m)
            y = pow(1. - x, self._m)
            dkdse = (1.0 - y) * ( (1. - y) / 2. + 2 * x * y / (1.-x)) / np.sqrt(se)
            return dkdse * dsdp / (1 - self._sr)
        elif pc <= 0:
            return 0.
        else:
            return 2 * pc * self._a + 3 * self._b * (pc**2)
