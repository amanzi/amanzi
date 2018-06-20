import numpy as np
#from scipy import optimize

import wrm_vangenuchten
import capillary_pressure
import permafrost_model_explicit_fpd

class PermafrostModel(object):
    def __init__(self, wrm=None, coef=1, dT=1.):
        if wrm is None:
            wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
        self.smooth_coef = coef
        self.dp = capillary_pressure.pc_ice(273.15-dT, 101325., 1)
        print self.dp
        self.pm = permafrost_model_explicit_fpd.PermafrostModel(wrm)

    def sm_coef(self, pcice, pcliq):
        mypcliq = max(pcliq, 0)
        dp = (pcice - mypcliq)/self.dp
        if 0 < dp < 1:
            return 1 + self.smooth_coef * (1-dp)*(1-dp)
        else:
            return 1
        
    def saturations_pc(self, pcliq, pcice):
        return self.pm.saturations_pc2(pcliq, pcice / self.sm_coef(pcice,pcliq), pcice)

    def saturations_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p,1)
        return self.saturations_pc(pcliq,pcice)
    
