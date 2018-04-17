import numpy as np

class BrooksCorey(object):
    def __init__( self, lambd=0., alpha=0., sr=0.0, smoothing_interval=0. ):
        self._lambda = lambd
        self._alpha = alpha
        self._sr = sr
        self._pc0 = smoothing_interval

        self._factor = -2.0 - (0.5 + 2.0) * self._lambda;
        self._pc_bubble = 1.0 / self._alpha
        
        if self._pc0 > 0.:
            k0 = self.k_relative(self._pc0) - 1.
            k0p = self.d_k_relative(self._pc0)
            self._a = (3 * k0 - k0p * self._pc0) / (self._pc0**2)
            self._b = (k0p * self._pc0 - 2 * k0) / (self._pc0**3)

    def capillaryPressure( self, s ):
        se = (s - self._sr) / (1.0 - self._sr)
        se = min(se, 1.0);
        return pow(se, -1.0/self._lambda) / self._alpha

    def d_capillaryPressure( self, s ):
        se = (s - self._sr) / (1.0 - self._sr)
        se = min(se, 1.0);
        return -1. / self._lambda * pow(se, -1.0/self._lambda - 1.) / self._alpha / (1. - self._sr);

    def saturation( self, pc ):
        if pc > self._pc_bubble:
            return pow(self._alpha * pc, -self._lambda) * (1.0 - self._sr) + self._sr
        else:
            return 1.0;

    def d_saturation( self, pc ):
        if pc > self._pc_bubble:
            return -pow(self._alpha * pc, -self._lambda - 1.0) * (1.0 - self._sr) * self._alpha * self._lambda
        else:
            return 0.

    def k_relative( self, pc ):
        if pc <= self._pc_bubble:
            return 1.0
        elif pc >= self._pc0:
            return pow(self._alpha * pc, self._factor)
        else:
            dpc = pc - self._pc_bubble
            return 1.0 + self._a * dpc**2 + self._b * dpc**3

    def d_k_relative( self, pc ):
        if pc <= self._pc_bubble:
            return 0.
        elif pc >= self._pc0:
            return self._factor * self._alpha * pow(self._alpha * pc, self._factor - 1.0)
        else:
            dpc = pc - self._pc_bubble
            return self._a * 2 * dpc + self._b * 3 * dpc**2
