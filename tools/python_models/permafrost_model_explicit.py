import numpy as np
#from scipy import optimize

import wrm_vangenuchten
import capillary_pressure

class PermafrostModel(object):
    def __init__(self, wrm=None):
        if wrm is None:
            wrm = wrm_vangenuchten.VanGenuchten(0.205, 5.5e-4)
        self.wrm = wrm

    #-----------

    def saturations_pc(self, pcliq, pcice):
        if pcliq <= 0.:
            # saturations if saturated
            sl = self.wrm.saturation(pcice)
            si = 1. - sl
            sg = 0.
            return sg, sl, si
        elif pcice <= pcliq:
            si = 0
            sl = self.wrm.saturation(pcliq)
            sg = 1 - sl
            return sg, sl, si
        else:
            sl = self.wrm.saturation(pcice)
            si = 1 - sl / self.wrm.saturation(pcliq)
            sg = 1 - sl - si
            return sg, sl, si

    def saturations_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        return self.saturations_pc(pcliq,pcice)

    #-----------

    def dsaturations_dpcliq_pc(self, pcliq, pcice):
        if pcliq <= 0.:
            # derivatives if saturated
            dsl = 0.
            dsi = 0.
            dsg = 0.
            return dsg, dsl, dsi

        elif pcice <= pcliq:
            dsi = 0.
            dsl = self.wrm.d_saturation(pcliq)
            dsg = -dsl
            return dsg, dsl, dsi
        else:
            dsl = 0.
            dsi = self.wrm.saturation(pcice) / (self.wrm.saturation(pcliq)**2) * self.wrm.d_saturation(pcliq)
            dsg = -dsl
            return dsg, dsl, dsi

    def dsaturations_dpcliq_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        return self.dsaturations_dpcliq_pc(pcliq, pcice)

    def dsaturations_dpcice_pc(self, pcliq, pcice):
        if pcliq <= 0.:
            # saturations if saturated
            dsl = self.wrm.d_saturation(pcice)
            dsi = -dsl
            dsg = 0.
            return dsg, dsl, dsi
        elif pcice <= pcliq:
            dsi = 0.
            dsl = 0.
            dsg = 0.
            return dsg, dsl, dsi
        else:
            dsl = self.wrm.d_saturation(pcice)
            dsi = -dsl / self.wrm.saturation(pcliq)
            dsg = - dsl - dsi
            return dsg, dsl, dsi

    def dsaturations_dpcice_Tp(self, T, p):
        pcliq = capillary_pressure.pc_liq(T,p)
        pcice = capillary_pressure.pc_ice(T,p)
        return self.dsaturations_dpcice_pc(pcliq, pcice)
