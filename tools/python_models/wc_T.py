import numpy as np
from scipy import optimize

import wrm_vangenuchten
import capillary_pressure
import eos

class WC_T(object):
    def __init__(self, pm):
        self.pm = pm

    def wc(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        sg,sl,si = self.pm.saturations_Tp(T,p)

        ni = eos.n_ice(T,p)
        nl = eos.n_water(T,p)
        ng = eos.n_gas(T,p)
        om = eos.mol_frac_gas(T)

        return sg*om*ng + sl*nl + si*ni

    def _wc_functor(self, T, wc):
        def functor(p):
            return wc - self.wc(T,p)
        return functor

    def pressure(self, T,wc):
        func = self._wc_functor(T,wc)

        dp = 101325
        a = 101325
        fa = func(a)
        if (fa < 0):
            fb = fa
            b = a
            nits = 0
            while fb <= 0 and nits < 1000:
                a = b
                b = b - dp
                fb = func(b)
                nits = nits + 1

            p,r = optimize.brentq(func,a,b, full_output=True)
            assert r.converged
        else:
            fb = fa
            b = a
            nits = 0
            while fb >= 0 and nits < 1000:
                a = b
                b = b + dp
                fb = func(b)
                nits = nits + 1

            p,r = optimize.brentq(func,b,a,full_output=True)
            assert r.converged

        return p
