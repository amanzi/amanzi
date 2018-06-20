import numpy as np
#from scipy import optimize

import wrm_vangenuchten
import capillary_pressure

class PermafrostModel(object):
    def __init__(self, wrm=None,
                 pc_ice=capillary_pressure.pc_ice, pc_liq=capillary_pressure.pc_liq):
        self.pc_ice = pc_ice
        self.pc_liq = pc_liq
        if wrm is None:
            wrm = wrm_vangenuchten.VanGenuchten(0.8, 1.5e-4)
        self.wrm = wrm

    #-----------

    def saturations_pc(self, pcliq, pcice):
        if pcliq <= 0.:
            # saturated
            sl = self.wrm.saturation(pcice)
            si = 1. - sl
            sg = 0.

        elif pcice <= pcliq:
            # unfrozen
            si = 0.
            sl = self.wrm.saturation(pcliq)
            sg = 1. - sl

        else:
            sl = self.wrm.saturation(pcice)
            si = 1 - sl / self.wrm.saturation(pcliq)
            sg = 1 - sl - si

        assert 0. <= sg <= 1.
        assert 0. <= sl <= 1.
        assert 0. <= si <= 1.
        return sg, sl, si

        
    def saturations_Tp(self, T, p):
        pcliq = self.pc_liq(T,p)
        pcice = self.pc_ice(T,p)
        return self.saturations_pc(pcliq,pcice)
